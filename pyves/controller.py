import json
import os
import sys
import _pyves
import numpy as np
import pandas as pd
import signal

from .signal_handler import SignalHandler



class Controller(object):
    def __init__(self) -> None:
        self.system = _pyves.System()
        self.signal_handler = SignalHandler()
        return



    def readForceField(self, path):
        with open(path, 'r') as ff_file:
            self.forcefield = json.loads(ff_file.read())
    


    def setForceField(self, ff:dict):
        self.forcefield = ff



    def readParameters(self, path):
        with open(path, 'r') as prms_file:
            _prms = json.loads(prms_file.read())
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
            # self.out_dir = _prms["output"]["dir"]
            # self.out_filename = _prms["output"]["filename"]
            # self.out_mode =_prms["output"]["mode"]

            # self.in_dir = _prms["input"]["dir"]
            # self.in_filename =_prms["input"]["filename"]
            # self.in_key = _prms["input"]["key"]



    def setParameters(self, prms:dict):
        self.parameters = prms
    


    def prepareSimulation(self):
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

        # Setup cells
        cells_per_dim = np.array(box_dims/self.cell_min_size).astype(int)
        cell_actual_size = box_dims/cells_per_dim
        for x in np.arange(0, box_dims[0], cell_actual_size[0]):
            for y in np.arange(0, box_dims[1], cell_actual_size[1]):
                for z in np.arange(0, box_dims[2], cell_actual_size[2]):
                    _min = np.array([x,y,z])
                    _max = np.array([x,y,z]) + cell_actual_size
                    self.system.cells.append(_pyves.Cell(_min, _max))
        for ci in self.system.cells:
            for cj in self.system.cells:
                if ci == cj:
                    ci.region.append(cj)
                if ci.isNeighbourOf(cj, self.system.box):
                    ci.proximity.append(cj)
                    ci.region.append(cj)
        
        # place particles in cells
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
        

        
    def sample(self, steps=None):
        if isinstance(steps, int):
            for i in range(steps):
                self.system.assertIntegrity()
                self.system.singleSimulationStep()
            self.time_actual += steps
            self.writeTrajectoryHDF() 
        else:
            self.time_actual = 0
            sampling_time_points = np.arange(self.time_actual, self.time_max, self.time_delta)
            for time_point in sampling_time_points:
                assert(self.time_actual == time_point)
                self.system.assertIntegrity()
                self.system.multipleSimulationSteps(self.time_delta)
                self.time_actual += self.time_delta
                self.writeTrajectoryHDF()



    def writeTrajectoryHDF(self, log=False):
        df = pd.DataFrame()
        for particle in self.system.particles:
            pos = self.system.box.scaleToBox(particle.position)
            data = dict(
                name = particle.name,
                x = pos[0],
                y = pos[1], 
                z = pos[2],
                ux = particle.ux,
                uy = particle.uy,
                uz = particle.uz,
                sigma = particle.sigma,
                epsilon = particle.epsilon,
                kappa = particle.kappa,
                gamma = particle.gamma
            )
            df = df.append(data, ignore_index=True, sort=False)
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
        # print(f"/time{self.time_actual} written")
        