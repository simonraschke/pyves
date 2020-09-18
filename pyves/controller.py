# pyves - python3 bindings for an easy use of vesicle2
# Copyright (C) 2020 Simon Raschke

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



import json
import os
from os import initgroups
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
from .point_distributions import *
from .analysis import analyze_df



class Controller(object):
    def __init__(self) -> None:
        SignalHandler.ProgramState = ProgramState.STARTUP
        self.system = _pyves.System()
        self.signal_handler = SignalHandler()
        self.time_actual = 0
        return



    # def readForceField(self, path):
    #     SignalHandler.ProgramState = ProgramState.SETUP
    #     with open(path, 'r') as ff_file:
    #         self.forcefield = json.loads(ff_file.read())
    


    # def setForceField(self, ff:dict):
    #     SignalHandler.ProgramState = ProgramState.SETUP
    #     self.forcefield = ff



    def readParameters(self, path):
        SignalHandler.ProgramState = ProgramState.SETUP
        with open(path, 'r') as prms_file:
            _prms = json.loads(prms_file.read())
            self.setParameters(_prms)



    def setParameters(self, _prms:dict):
        if not SignalHandler.ProgramState == ProgramState.SHUTDOWN:
            SignalHandler.ProgramState = ProgramState.SETUP
        self.prms_complete = _prms
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

        self.system.interaction_cutoff = _prms["system"]["interaction"]["cutoff"]

        self.system.update_interval = _prms["control"]["update_interval"]
        self.system.neighbor_cutoff = _prms["control"]["neighbor_cutoff"]

        self.time_delta = _prms["control"]["time_delta"]
        self.time_max = _prms["control"]["time_max"]
        self.particle_prms = _prms["system"]["particles"]
        self.cell_min_size = _prms["control"]["cell_min_size"]

        self.output = _prms["output"]
        self.input = _prms["input"]

        assert(self.system.neighbor_cutoff >= self.system.interaction_cutoff)
        assert(self.cell_min_size >= self.system.interaction_cutoff)
    


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

    

    def makeInteractionLookupTable(self):
        particles = _pyves.ParticleContainer()
        for name, particle_ff in self.particle_prms.items():
            particles.append(_pyves.Particle([0,0,0], [0,1,0], 
                sigma=particle_ff["sigma"], 
                kappa=particle_ff["kappa"], 
                eps=particle_ff["epsilon"], 
                gamma=particle_ff["gamma"], 
                name=name)
            )
        self.system.makeLookupTableFrom(particles)
        import pprint
        pprint.pprint(self.system.lookupTable)



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
        for name, particle_ff in self.particle_prms.items():
            
            if "sphere" in particle_ff["dist"]:
                r0 = 2.0**(1.0/6)*particle_ff["sigma"]
                points = sunflower_sphere_points(particle_ff["number"])
                if particle_ff["dist"] == "sphere":
                    radius = r0/(2.0*np.sin(particle_ff["gamma"]))
                else:
                    optimum_size = int(re.findall('(?<=sphere)\d+', particle_ff["dist"])[0])
                    radius = r0/4*np.sqrt(1.1027*optimum_size)
                for point in points:
                    self.system.particles.append(_pyves.Particle(point*radius+box_dims/2, point, 
                        sigma=particle_ff["sigma"], kappa=particle_ff["kappa"], eps=particle_ff["epsilon"], gamma=particle_ff["gamma"], name=name))
                    if particle_ff["bound_translation"] != None:
                        self.system.particles[-1].translation_bound_sq  = particle_ff["bound_translation"]**2
                    if particle_ff["bound_rotation"] != None:
                        self.system.particles[-1].rotation_bound = particle_ff["bound_rotation"]
            
            elif "plane" in particle_ff["dist"]:
                points = grid_plane_points(particle_ff["number"])
                area = float(re.findall('(?<=plane)\d+', particle_ff["dist"])[0])
                edge = np.sqrt(area)
                for point in points:
                    self.system.particles.append(_pyves.Particle(point*edge/2+box_dims/2, [0,0,1], 
                        sigma=particle_ff["sigma"], kappa=particle_ff["kappa"], eps=particle_ff["epsilon"], gamma=particle_ff["gamma"], name=name))
                    if particle_ff["bound_translation"] != None:
                        self.system.particles[-1].translation_bound_sq  = particle_ff["bound_translation"]**2
                    if particle_ff["bound_rotation"] != None:
                        self.system.particles[-1].rotation_bound = particle_ff["bound_rotation"]

            elif particle_ff["dist"] == "random":
                for _ in range(particle_ff["number"]):
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
                    if particle_ff["bound_translation"] != None:
                        self.system.particles[-1].translation_bound_sq  = particle_ff["bound_translation"]**2
                    if particle_ff["bound_rotation"] != None:
                        self.system.particles[-1].rotation_bound = particle_ff["bound_rotation"]
            else:
                raise NotImplementedError(f"Particle distribution {particle_ff['dist']} not implemented.")

        self.setupCells()
        self.placeParticlesInCells()
        self.writeTrajectoryHDF(timestats=True)



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
            self.system.particles[-1].initial_position = [row["initial_x"], row["initial_y"], row["initial_z"]]
            self.system.particles[-1].initial_orientation = [row["initial_ux"], row["initial_uy"], row["initial_uz"]]
            self.system.particles[-1].translation_bound_sq = row["translation_bound_sq"]
            self.system.particles[-1].rotation_bound = row["rotation_bound"]
        
        self.setupCells()
        self.placeParticlesInCells()



    def prepareSimulation(self):
        if not SignalHandler.ProgramState == ProgramState.SHUTDOWN:
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
            self.makeInteractionLookupTable()
            self.system.makeNeighborLists()



    def sample(self, steps=None, timestats=False):
        if not SignalHandler.ProgramState == ProgramState.SHUTDOWN:
            SignalHandler.ProgramState = ProgramState.RUNNING

            self.sample_starttime = time.perf_counter()

            if isinstance(steps, int):
                for i in range(steps):
                    if not SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                        self.system.assertIntegrity()
                        self.system.singleSimulationStep()
                        self.time_actual += 1
                    else:
                        break
                if timestats:
                    self.sample_endtime = time.perf_counter()
                    print(
                        f"time {self.time_actual} took {self.sample_endtime-self.sample_starttime:.4f} s", " | ", 
                        f"{(self.sample_endtime-self.sample_starttime)*1000*1000/len(self.system.particles)/steps:.4f} ns /particle/step", " | ", 
                        end=""
                    )
                # print("writing...")
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
                            self.sample_endtime = time.perf_counter()
                            print(
                                f"sampling to step {self.time_actual} took {self.sample_endtime-self.sample_starttime:.4f} s", " | ", 
                                f"{(self.sample_endtime-self.sample_starttime)*1000*1000/len(self.system.particles)/self.time_delta:.4f} ns /particle/step", " | ", 
                                end=""
                            )
                        self.writeTrajectoryHDF(timestats=timestats)
                        if timestats: print()

                    if SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                        print("shutting down")
                        sys.exit(0)
                    self.sample_starttime = time.perf_counter()



    def writeTrajectoryHDF(self, log=False, timestats=False):
        starttime = time.perf_counter()

        positions = np.array([self.system.box.scaleToBox(p.position) for p in self.system.particles])
        initial_positions = np.array([self.system.box.scaleToBox(p.initial_position) for p in self.system.particles])
        initial_orientations = np.array([self.system.box.scaleToBox(p.initial_orientation) for p in self.system.particles])

        df = pd.DataFrame(dict(
            name = [p.name for p in self.system.particles],
            x = [p[0] for p in positions],
            y = [p[1] for p in positions],
            z = [p[2] for p in positions],
            ux = [p.ux for p in self.system.particles],
            uy = [p.uy for p in self.system.particles],
            uz = [p.uz for p in self.system.particles],
            initial_x = [p[0] for p in initial_positions],
            initial_y = [p[1] for p in initial_positions],
            initial_z = [p[2] for p in initial_positions],
            translation_bound_sq = [p.translation_bound_sq for p in self.system.particles],
            initial_ux = [p[0] for p in initial_orientations],
            initial_uy = [p[1] for p in initial_orientations],
            initial_uz = [p[2] for p in initial_orientations],
            rotation_bound = [p.rotation_bound for p in self.system.particles],
            sigma = [p.sigma for p in self.system.particles],
            epsilon = [p.epsilon for p in self.system.particles],
            kappa = [p.kappa for p in self.system.particles],
            gamma = [p.gamma for p in self.system.particles]
        ))

        df.to_hdf(
            path_or_buf=os.path.join(self.output["dir"], self.output["filename"]),
            key=f"/time{self.time_actual}",
            mode="a",# if self.output["mode"] == "append" else "w",
            # append=True if self.output["mode"] == "append" else False,
            format="table",
            complevel=1,
        )

        node = h5py.File(os.path.join(self.output["dir"], self.output["filename"]), mode="a").get(f"/time{self.time_actual}")
        node.attrs["box.x"] = self.system.box.x
        node.attrs["box.y"] = self.system.box.y
        node.attrs["box.z"] = self.system.box.z
        node.attrs["threads"] = self.system.threads
        node.attrs["temperature"] = self.system.temperature
        node.attrs["translation_step"] = self.system.translation()
        node.attrs["rotation_step"] = self.system.rotation()
        node.attrs["interaction.cutoff"] = self.system.interaction_cutoff
        node.attrs["time"] = self.time_actual
        node.attrs["cell_min_size"] = self.prms_complete["control"]["cell_min_size"]
        try:
            node.attrs["simulation_walltime"] = time.perf_counter() - self.sample_starttime
        except:
            node.attrs["simulation_walltime"] = time.perf_counter()

        if log:
            print(df)
            print(df.info())
            
        if timestats:
            endtime = time.perf_counter()
            print(f" write HDF5 took {endtime-starttime:.4f} s", end="")
        