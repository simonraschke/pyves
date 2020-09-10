import json
import os
from sqlite3.dbapi2 import TimestampFromTicks
import sys
import _pyves
import numpy as np
import pandas as pd
import signal
import time
import h5py
import re

from .signal_handler import SignalHandler, ProgramState



class Controller(object):
    def __init__(self) -> None:
        SignalHandler.ProgramState = ProgramState.STARTUP
        self.system = _pyves.System()
        self.signal_handler = SignalHandler()
        self.time_actual = 0
        return



    def readForceField(self, path):
        SignalHandler.ProgramState = ProgramState.SETUP
        with open(path, 'r') as ff_file:
            self.forcefield = json.loads(ff_file.read())
    


    def setForceField(self, ff:dict):
        SignalHandler.ProgramState = ProgramState.SETUP
        self.forcefield = ff



    def readParameters(self, path):
        SignalHandler.ProgramState = ProgramState.SETUP
        with open(path, 'r') as prms_file:
            _prms = json.loads(prms_file.read())
            self.setParameters(_prms)



    def setParameters(self, _prms:dict):
        SignalHandler.ProgramState = ProgramState.SETUP
        self.system.cores = _prms["hardware"]["cores"]
        self.system.threads = _prms["hardware"]["threads"]

        self.system.temperature = _prms["system"]["temperature"]
        self.system.box.x = _prms["system"]["box"]["x"]
        self.system.box.y = _prms["system"]["box"]["y"]
        self.system.box.z = _prms["system"]["box"]["z"]
        translation = _prms["system"]["translation"]
        self.system.translation = _pyves.StepwidthAlignmentUnit()
        self.system.translation.setup(translation["interval"], translation["min"], translation["max"], translation["target"])
        rotation = _prms["system"]["rotation"]
        self.system.rotation = _pyves.StepwidthAlignmentUnit()
        self.system.rotation.setup(rotation["interval"], rotation["min"], rotation["max"], rotation["target"])

        self.time_delta = _prms["control"]["time_delta"]
        self.time_max = _prms["control"]["time_max"]
        self.particle_prms = _prms["control"]["particles"]
        self.cell_min_size = _prms["control"]["cell_min_size"]

        self.output = _prms["output"]
        self.input = _prms["input"]
    


    def setupCells(self):
        # Setup cells
        box_dims = np.array([self.system.box.x, self.system.box.y, self.system.box.z])
        cells_per_dim = np.array(box_dims/self.cell_min_size).astype(int)
        cell_actual_size = box_dims/cells_per_dim

        for x in np.arange(0, box_dims[0], cell_actual_size[0]):
            for y in np.arange(0, box_dims[1], cell_actual_size[1]):
                for z in np.arange(0, box_dims[2], cell_actual_size[2]):
                    _min = np.array([x,y,z])
                    _max = np.array([x,y,z]) + cell_actual_size
                    self.system.cells.append(_pyves.Cell(_min, _max, box=self.system.box))

        for ci in self.system.cells:
            for cj in self.system.cells:
                if ci == cj:
                    ci.region.append(cj)
                if ci.isNeighbourOf(cj):
                    ci.proximity.append(cj)
                    ci.region.append(cj)
        
        for c in self.system.cells:
            assert(c.assertIntegrity())

    

    def placeParticlesInCells(self):
        cell_place_counter = 0
        for i, particle in enumerate(self.system.particles):
            for j, cell in enumerate(self.system.cells):
                if cell.insideCellBounds(particle):
                    cell.particles.append(particle)
                    # print(i, "to", j)
                    cell_place_counter += 1
                    break
        assert(cell_place_counter == len(self.system.particles))
        assert(self.system.assertIntegrity())



    def prepareNew(self):
        box_dims = np.array([self.system.box.x, self.system.box.y, self.system.box.z])

        # distribute particles
        for name, particle_prms in self.particle_prms.items():
            assert( name in self.forcefield )
            particle_ff = self.forcefield[name]
            if particle_prms["dist"] != "random":
                raise NotImplementedError("no fixed distribution supported")
            elif particle_prms["dist"] == "random":
                for _ in range(particle_prms["number"]):
                    self.system.particles.append(_pyves.Particle(np.random.rand(3)*box_dims, np.random.uniform(-1,1,3), 
                        sigma=particle_ff["sigma"], kappa=particle_ff["kappa"], eps=particle_ff["epsilon"], gamma=particle_ff["gamma"], name=name))
                    # repeat until particle is free
                    particle_try_set_count = 0
                    while not self.system.particleIsFree(self.system.particles[-1]):
                        particle_try_set_count += 1
                        self.system.particles[-1].position = np.random.rand(3)*box_dims
                        if particle_try_set_count > 1e6:
                            print(f"tried to set particle {particle_try_set_count} time. Aborting")
                            sys.exit(signal.SIGKILL)

        self.setupCells()
        self.placeParticlesInCells()



    def prepareFromData(self, path):
        if self.input["key"].lower() == "head":
            with h5py.File(path, 'r') as h5file:
                key = sorted([s for s in h5file.keys() if s.startswith("time")], key=lambda x:int(re.findall('(?<=time)\d+', x)[0]))[-1]
        else:
            key = self.input["key"]

        df = pd.read_hdf(path_or_buf=path, key=key)
        
        self.time_actual = int(re.findall('(?<=time)\d+', key)[0])

        for _, row in df.iterrows():
            self.system.particles.append(_pyves.Particle([row["x"], row["y"], row["z"]], [row["ux"], row["uy"], row["uz"]], 
                sigma=row["sigma"], kappa=row["kappa"], eps=row["epsilon"], gamma=row["gamma"], name=row["name"]))
        
        self.setupCells()
        self.placeParticlesInCells()



    def prepareSimulation(self):
        SignalHandler.ProgramState = ProgramState.SETUP

        try:
            input_data_file_path = os.path.join(self.input["dir"], self.input["filename"])
            self.prepareFromData(input_data_file_path)
        except TypeError as e:
            print("TypeError: unable to read input file:", e)
            print("creating System new")
            self.prepareNew()
        except OSError as e:
            print("TypeError: unable to read datafile:", e)
            print("creating System new")
            self.prepareNew()



    def sample(self, steps=None, timestats=False):
        SignalHandler.ProgramState = ProgramState.RUNNING
        starttime = time.perf_counter()

        if isinstance(steps, int):
            for i in range(steps):
                if not SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                    self.system.assertIntegrity()
                    self.system.singleSimulationStep()
                    self.time_actual += 1
                else:
                    break
            if timestats:
                endtime = time.perf_counter()
                print(
                    f"time {self.time_actual} took {endtime-starttime:.4f} s", " | ", 
                    f"{(endtime-starttime)*1000*1000/len(self.system.particles)/steps:.4f} ns /particle/step", " | ", 
                    end=""
                )
            self.writeTrajectoryHDF(timestats=timestats)
            if timestats: print()
            if SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                print("shutting down")
                sys.exit(0)

        else:
            sampling_time_points = np.arange(self.time_actual, self.time_max, self.time_delta)
            for time_point in sampling_time_points:
                if not SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                    assert(self.time_actual == time_point)
                    self.system.assertIntegrity()
                    self.system.multipleSimulationSteps(self.time_delta)
                    self.time_actual += self.time_delta
                    if timestats:
                        endtime = time.perf_counter()
                        print(
                            f"sampling to step {self.time_actual} took {endtime-starttime:.4f} s", " | ", 
                            f"{(endtime-starttime)*1000*1000/len(self.system.particles)/self.time_delta:.4f} ns /particle/step", " | ", 
                            end=""
                        )
                    self.writeTrajectoryHDF(timestats=timestats)
                    if timestats: print()

                if SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                    print("shutting down")
                    sys.exit(0)
                starttime = time.perf_counter()



    def writeTrajectoryHDF(self, log=False, timestats=False):
        starttime = time.perf_counter()

        positions = np.array([self.system.box.scaleToBox(p.position) for p in self.system.particles])

        df = pd.DataFrame(dict(
            name = [p.name for p in self.system.particles],
            x = [p[0] for p in positions],
            y = [p[1] for p in positions],
            z = [p[2] for p in positions],
            ux = [p.ux for p in self.system.particles],
            uy = [p.uy for p in self.system.particles],
            uz = [p.uz for p in self.system.particles],
            sigma = [p.sigma for p in self.system.particles],
            epsilon = [p.epsilon for p in self.system.particles],
            kappa = [p.kappa for p in self.system.particles],
            gamma = [p.gamma for p in self.system.particles]
        ))

        df.to_hdf(
            path_or_buf=os.path.join(self.output["dir"], self.output["filename"]),
            key=f"/time{self.time_actual}",
            mode="r+" if self.output["mode"] == "append" else "a",
            append=True if self.output["mode"] == "append" else False,
            format="table",
            complevel=1,
        )

        if log:
            print(df)
            print(df.info())
            
        if timestats:
            endtime = time.perf_counter()
            print(f" write HDF5 took {endtime-starttime:.4f} s", end="")
        