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
import sys
import _pyves
import numpy as np
import pandas as pd
import signal
import time
import h5py
import re
import itertools

from more_itertools import windowed

from ._version import __version__
from .signal_handler import SignalHandler, ProgramState
from .point_distributions import *
from .analysis import analyzeSnapshot, analyzeTrajectory
from .utility import h5store, h5load



class Controller():
    def __init__(self) -> None:
        SignalHandler.ProgramState = ProgramState.STARTUP
        self.system = _pyves.System()
        self.signal_handler = SignalHandler()
        self.time_actual = 0
        return


    
    @staticmethod
    def printRuntimeInfo():
        print("pyves version:", __version__)
        print("pyves concurrency model:", _pyves.concurrency_model())
        condaenv = os.environ.get("CONDA_DEFAULT_ENV", None)
        if not isinstance(condaenv, type(None)):
            print("conda environment:", condaenv)
        virtualenv = os.environ.get("VIRTUAL_ENV", None)
        if not isinstance(virtualenv, type(None)):
            print("virtual environment:", virtualenv)
        if isinstance(condaenv, type(None)) and isinstance(virtualenv, type(None)):
            print("running in default python environment:", sys.executable)



    @classmethod
    def Static(
        cls, 
        prmspath,
        timestats = True,
        analysis = True,
        analysis_inline = False,
        benchmark = True
    ):
        Controller.printRuntimeInfo()
        ctrl = cls()
        ctrl.readParameters(prmspath)
        ctrl.prepareSimulation()

        if benchmark:
            ctrl.system.benchmark(100)

        ctrl.sample(timestats=timestats, analysis=analysis_inline)

        if analysis and not analysis_inline:
            analyzeTrajectory(prmspath=prmspath, timestats=timestats, threads=-1)
        return ctrl



    @classmethod
    def DynamicSystemFlow(
        cls, 
        prmspath,
        attr,
        times,
        values,
        timestats = True,
        analysis = True,
        analysis_inline = False,
        benchmark = True
    ):
        assert len(times) == np.unique(times).size

        Controller.printRuntimeInfo()
        ctrl = cls()
        ctrl.readParameters(prmspath)
        # ctrl.time_max = max(times)
        ctrl.prepareSimulation()

        if benchmark:
            ctrl.system.benchmark(100)
        
        assert hasattr(ctrl.system, attr), attr

        times = np.array(times, dtype=np.int64)
        values = np.array(values)
        time_indices = times.argsort()
        times = times[time_indices]
        values = values[time_indices]

        print(*zip(times, values))

        if ctrl.time_actual < min(times):
            ctrl.sample(steps=min(times)-ctrl.time_actual, timestats=timestats, analysis=analysis_inline)
        setattr(ctrl.system, attr, values[times.tolist().index(ctrl.time_actual)])

        for start, end in windowed(times, 2):
            if ctrl.time_actual < end and ctrl.time_actual >= start:
                ctrl.sample(steps=end-ctrl.time_actual, timestats=timestats, analysis=analysis_inline)
                assert ctrl.time_actual == end, ctrl.time_actual
                assert end in times, end
                assert hasattr(ctrl.system, attr), attr

                setattr(ctrl.system, attr, values[times.tolist().index(ctrl.time_actual)])

                assert abs(getattr(ctrl.system, attr) - values[times.tolist().index(ctrl.time_actual)]) < 1e-5, abs(getattr(ctrl.system, attr) - values[times.tolist().index(ctrl.time_actual)])
                assert abs(ctrl.getMetadata()[attr] - values[times.tolist().index(ctrl.time_actual)]) < 1e-5, abs(ctrl.getMetadata()[attr] - values[times.tolist().index(ctrl.time_actual)])                
                assert attr in h5load(os.path.join(ctrl.output["dir"], ctrl.output["filename"]), f"time{ctrl.time_actual}")[1]
                assert abs(h5load(os.path.join(ctrl.output["dir"], ctrl.output["filename"]), f"time{ctrl.time_actual}")[1][attr] - values[times.tolist().index(ctrl.time_actual)-1]) < 1e-5

        if analysis and not analysis_inline:
            analyzeTrajectory(prmspath=prmspath, timestats=timestats, threads=-1)

        return ctrl



    @classmethod
    def DynamicParticleFlow(
        cls, 
        prmspath,
        attr,
        times,
        values,
        timestats = True,
        analysis = True,
        analysis_inline = False,
    ):
        assert len(times) == np.unique(times).size

        Controller.printRuntimeInfo()
        ctrl = cls()
        ctrl.readParameters(prmspath)
        # ctrl.time_max = max(times)
        ctrl.prepareSimulation()
        
        for particle in ctrl.system.particles:
            assert hasattr(particle, attr), attr

        times = np.array(times, dtype=np.int64)
        values = np.array(values)
        time_indices = times.argsort()
        times = times[time_indices]
        values = values[time_indices]

        print(*zip(times, values))

        if ctrl.time_actual < min(times):
            ctrl.sample(steps=min(times)-ctrl.time_actual, timestats=timestats, analysis=analysis_inline)

        # for particle in ctrl.system.particles:
        #     setattr(particle, attr, values[times.tolist().index(ctrl.time_actual)])
        for i, _ in enumerate(ctrl.system.particles):
        #     print("HHHHHHHIIIIIIIIIIIIIEEEEEEEEEEEEEERRRRRRRRRR", i)
        #     print("HHHHHHHIIIIIIIIIIIIIEEEEEEEEEEEEEERRRRRRRRRR", ctrl.system.particles[i])
        #     print("HHHHHHHIIIIIIIIIIIIIEEEEEEEEEEEEEERRRRRRRRRR", values)
        #     print("HHHHHHHIIIIIIIIIIIIIEEEEEEEEEEEEEERRRRRRRRRR", times)
        #     print("HHHHHHHIIIIIIIIIIIIIEEEEEEEEEEEEEERRRRRRRRRR", times.tolist())
        #     print("HHHHHHHIIIIIIIIIIIIIEEEEEEEEEEEEEERRRRRRRRRR", ctrl.time_actual)
        #     print("HHHHHHHIIIIIIIIIIIIIEEEEEEEEEEEEEERRRRRRRRRR", times.tolist().index(ctrl.time_actual))
            setattr(ctrl.system.particles[i], attr, values[times.tolist().index(ctrl.time_actual)])

        for start, end in windowed(times, 2):
            if ctrl.time_actual < end and ctrl.time_actual >= start:
                ctrl.sample(steps=end-ctrl.time_actual, timestats=timestats, analysis=analysis_inline)
                assert ctrl.time_actual == end, ctrl.time_actual
                assert end in times, end
                for particle in ctrl.system.particles:
                    assert hasattr(particle, attr), attr

                for i, _ in enumerate(ctrl.system.particles):
                    setattr(ctrl.system.particles[i], attr, values[times.tolist().index(ctrl.time_actual)])

                _temp_traj_data = h5load(os.path.join(ctrl.output["dir"], ctrl.output["filename"]), f"time{ctrl.time_actual}")[0]
                for i, particle in enumerate(ctrl.system.particles):
                    assert abs(getattr(particle, attr) - values[times.tolist().index(ctrl.time_actual)]) < 1e-5, abs(getattr(particle, attr) - values[times.tolist().index(ctrl.time_actual)])
                    print(getattr(particle, attr), _temp_traj_data.iloc[i][attr])
                    assert abs(getattr(particle, attr) - _temp_traj_data.iloc[i][attr]) < 1e-5
                    
                # assert abs(h5load(os.path.join(ctrl.output["dir"], ctrl.output["filename"]), f"time{ctrl.time_actual}")[1][attr] - values[times.tolist().index(ctrl.time_actual)-1]) < 1e-5

        if analysis and not analysis_inline:
            analyzeTrajectory(prmspath=prmspath, timestats=timestats, threads=-1)

        return ctrl



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

        self.system.exchange_global_number = _prms["system"]["exchange"]["global"].get("number", 0)
        self.system.exchange_global_etot_theshold = _prms["system"]["exchange"]["global"].get("etot_threshold", -1.0)
        self.system.exchange_global_orientation = _prms["system"]["exchange"]["global"].get("orientation", False)
        self.system.exchange_local_number = _prms["system"]["exchange"]["local"].get("number", 0)
        self.system.exchange_local_etot_theshold = _prms["system"]["exchange"]["local"].get("etot_threshold", -1.0)
        self.system.exchange_local_orientation = _prms["system"]["exchange"]["local"].get("orientation", False)        

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
        self.system.interaction_surface = _prms["system"]["interaction"].get("z_surface", False)
        self.system.interaction_surface_width = _prms["system"]["interaction"].get("z_surface_width", 0)

        self.system.cell_update_interval = _prms["control"]["cell_update_interval"]
        self.system.neighbor_update_interval = _prms["control"]["neighbor_update_interval"]
        self.system.neighbor_cutoff = _prms["control"]["neighbor_cutoff"]

        self.time_delta = _prms["control"]["time_delta"]
        self.time_max = _prms["control"]["time_max"]
        self.particle_prms = _prms["system"]["particles"]
        self.cell_min_size = _prms["control"]["cell_min_size"]

        self.output = _prms["output"]
        self.input = _prms["input"]

        self.direct_analysis = self.output.get("direct_analysis", False)

        assert(self.system.neighbor_cutoff >= self.system.interaction_cutoff)
        assert(self.cell_min_size >= self.system.interaction_cutoff)
    


    def setForAllParticles(self, attr=dict(), condition=dict()):
        assert isinstance(attr, dict)
        assert isinstance(condition, dict)
        print("setting all", condition, "particles to", attr)
        for particle in self.system.particles:
            if sum([ getattr(particle, k) == v for k,v in condition.items() ]) == len(condition.items()):
                for k, v in attr.items():
                    setattr(particle, k, v)



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
        
        print(f"generated {len(self.system.cells)} cells")



    # def makeInteractionLookupTable(self, force_zero=False):
    #     particles = _pyves.ParticleContainer()
    #     for name, particle_ff in self.particle_prms.items():
    #         if not force_zero and particle_ff.get("number", particle_ff["ratio"]) == 0:
    #             continue
    #         particles.append(
    #             _pyves.Particle(
    #                 [0,0,0], [0,1,0], 
    #                 sigma=particle_ff["sigma"], 
    #                 kappa=particle_ff["kappa"], 
    #                 eps=particle_ff["epsilon"], 
    #                 gamma=particle_ff["gamma"], 
    #                 name=name
    #             )
    #         )
    #     self.system.makeLookupTableFrom(particles)
    #     # import pprint
    #     # pprint.pprint(self.system.lookupTable)



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

        data = [len(cell.particles) for cell in self.system.cells]
        unique, counts = np.unique(data, return_counts=True)
        for u, c in zip(unique, counts):
            print(f"cells with {u:>2} particles: {c:>5}")



    def prepareNew(self):
        def _additional_particle_setup(particle, ff):
            if self.system.interaction_surface:
                particle.surface_affinity_translation = ff["surface_affinity_translation"]
                particle.surface_affinity_rotation = ff["surface_affinity_rotation"]
            if "self_affinity" in ff.keys():
                particle.self_affinity = ff["self_affinity"]
            if "other_affinity" in ff.keys():
                particle.other_affinity = ff["other_affinity"]
            if ff["bound_translation"] != None:
                particle.translation_bound_sq  = ff["bound_translation"]**2
            if ff["bound_rotation"] != None:
                particle.rotation_bound = ff["bound_rotation"]

        box_dims = np.array([self.system.box.x, self.system.box.y, self.system.box.z])
        p_add_counter = 0

        count_guv = sum([1 for _, particle_ff in self.particle_prms.items() if "guv" in particle_ff["dist"]])
        if count_guv > 0:
            tuples = np.array([[k,v["number"],v["sigma"]] for k,v in self.particle_prms.items() if "guv" in v["dist"]])
            possible_names = tuples[:,0]
            numbers = np.array(tuples[:,1], dtype="int")
            possibilities = numbers / np.sum(numbers)
            sigmas = np.array(tuples[:,2], dtype="float")

            r0s = sigmas*2.0**(1.0/6)
            areas_per_ptype = 2.0 * np.sqrt(3.0) * (r0s/2)**2*numbers
            area_guv = np.sum(areas_per_ptype)
            radius = np.sqrt(area_guv/(np.pi*4))
            points = sunflower_sphere_points(int(np.sum(numbers)))
            namearray = list(itertools.chain.from_iterable(itertools.repeat(name, number) for name, number in zip(possible_names, numbers)))
            print(f"setting up a GUV with the radius = {radius:.3f} from {int(np.sum(numbers))} particles")

            assert len(points) == np.sum(numbers)
            assert all(radius < box_dims/2)
            
            if sum([1 for _,v in self.particle_prms.items() if "guvseparated" in v["dist"]]):
                print("generating a sorted GUV")
                namearray = sorted(namearray)
            else:
                import random
                value_changes = []
                N = 500
                print(f"checking {N} possible combinations for largest variation")
                for _ in range(N):
                    random.shuffle(namearray)
                    value_changes.append((np.where(np.roll(namearray,1)!=namearray)[0].size, namearray))
                value_changes = np.array(value_changes, dtype='object')
                largest_var = np.argmax(value_changes[:,0].astype(int))
                print(f"found largest variation in {largest_var} with {value_changes[largest_var,0]} value changes")
                namearray = value_changes[np.argmax(value_changes[:,0].astype(int)),1]
                del value_changes

            for name, point in zip(namearray, points):
                particle_ff = self.particle_prms[name]
                self.system.particles.append(_pyves.Particle(point*radius+(box_dims/2), point, 
                    sigma=particle_ff["sigma"], kappa=particle_ff["kappa"], eps=particle_ff["epsilon"], gamma=particle_ff["gamma"], name=name))
                p_add_counter += 1
                # print(p_add_counter, "added", self.system.particles[-1])
                _additional_particle_setup(self.system.particles[-1], particle_ff)

        count_hexplane = sum([1 for _, particle_ff in self.particle_prms.items() if "hexplane" in particle_ff["dist"]])
        if count_hexplane > 0:
            average_sigma = sum([particle_ff["sigma"] for _, particle_ff in self.particle_prms.items() if "hexplane" in particle_ff["dist"]])/count_hexplane
            points = hexagonal_lattice_points(x=box_dims[0], y=box_dims[1], distance=average_sigma*2.0**(1.0/6))
            possible_names = [k for k,v in self.particle_prms.items() if "hexplane" in v["dist"]]
            possibilities = np.array([v["ratio"] for _,v in self.particle_prms.items() if "hexplane" in v["dist"]])
            possibilities = possibilities / np.cumsum(possibilities)[-1]
            namearray = np.random.choice(possible_names, len(points), p=possibilities)

            for name, point in zip(namearray, points):
                particle_ff = self.particle_prms[name]
                point[2] = box_dims[2]/2
                self.system.particles.append(_pyves.Particle(point, [0,0,1], 
                    sigma=particle_ff["sigma"], kappa=particle_ff["kappa"], eps=particle_ff["epsilon"], gamma=particle_ff["gamma"], name=name))
                p_add_counter += 1
                # print(p_add_counter, "added", self.system.particles[-1])
                _additional_particle_setup(self.system.particles[-1], particle_ff)
        
        # distribute particles
        for name, particle_ff in self.particle_prms.items():
            if "guv" in particle_ff["dist"]:
                continue

            elif "sphere" in particle_ff["dist"]:
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
                    p_add_counter += 1
                    # print(p_add_counter, "added", self.system.particles[-1])
                _additional_particle_setup(self.system.particles[-1], particle_ff)
                    # if self.system.interaction_surface:
                    #     self.system.particles[-1].surface_affinity_translation = particle_ff["surface_affinity_translation"]
                    #     self.system.particles[-1].surface_affinity_rotation = particle_ff["surface_affinity_rotation"]
                    # if particle_ff["bound_translation"] != None:
                    #     self.system.particles[-1].translation_bound_sq  = particle_ff["bound_translation"]**2
                    # if particle_ff["bound_rotation"] != None:
                    #     self.system.particles[-1].rotation_bound = particle_ff["bound_rotation"]

            
            elif "plane" in particle_ff["dist"]:
                if "hex" in particle_ff["dist"]:
                    continue
                else:
                    points = grid_plane_points(particle_ff["number"])
                    area = float(re.findall('(?<=plane)\d+', particle_ff["dist"])[0])
                    edge = np.sqrt(area)
                    for point in points:
                        self.system.particles.append(_pyves.Particle(point*edge/2+box_dims/2, [0,0,1], 
                            sigma=particle_ff["sigma"], kappa=particle_ff["kappa"], eps=particle_ff["epsilon"], gamma=particle_ff["gamma"], name=name))
                        p_add_counter += 1
                        # print(p_add_counter, "added", self.system.particles[-1])
                        _additional_particle_setup(self.system.particles[-1], particle_ff)
                            
                        # if self.system.interaction_surface:
                        #     self.system.particles[-1].surface_affinity_translation = particle_ff["surface_affinity_translation"]
                        #     self.system.particles[-1].surface_affinity_rotation = particle_ff["surface_affinity_rotation"]
                        # if particle_ff["bound_translation"] != None:
                        #     self.system.particles[-1].translation_bound_sq  = particle_ff["bound_translation"]**2
                        # if particle_ff["bound_rotation"] != None:
                        #     self.system.particles[-1].rotation_bound = particle_ff["bound_rotation"]

            elif particle_ff["dist"] == "random":
                for i in range(particle_ff["number"]):
                    self.system.particles.append(_pyves.Particle(np.random.rand(3)*box_dims, np.random.uniform(-1,1,3), 
                        sigma=particle_ff["sigma"], kappa=particle_ff["kappa"], eps=particle_ff["epsilon"], gamma=particle_ff["gamma"], name=name))
                    p_add_counter += 1
                    # print(p_add_counter, "added", self.system.particles[-1])
                    if self.system.interaction_surface:
                        self.system.particles[-1].surface_affinity_translation = particle_ff["surface_affinity_translation"]
                        self.system.particles[-1].surface_affinity_rotation = particle_ff["surface_affinity_rotation"]

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

        print(f"generated {len(self.system.particles)} particles")

        self.setupCells()
        self.placeParticlesInCells()
        self.writeTrajectoryHDF(timestats=True)



    def prepareFromData(self, path, df = None):
        if isinstance(df, type(None)):
            if self.input["key"].lower() == "head":
                with h5py.File(path, 'r') as h5file:
                    key = sorted([s for s in h5file.keys() if s.startswith("time")], key=lambda x:int(re.findall('(?<=time)\d+', x)[0]))[-1]
            else:
                key = self.input["key"]

            # df, metadata = pd.read_hdf(path_or_buf=path, key=key)
            df, metadata = h5load(path, key)
            self.time_actual = int(re.findall('(?<=time)\d+', key)[0])
            self.setFromMetadata(metadata)

        for _, row in df.iterrows():
            self.system.particles.append(_pyves.Particle([row["x"], row["y"], row["z"]], [row["ux"], row["uy"], row["uz"]], 
                sigma=row["sigma"], kappa=row["kappa"], eps=row["epsilon"], gamma=row["gamma"], name=row["name"]))
            self.system.particles[-1].initial_position = [row["initial_x"], row["initial_y"], row["initial_z"]]
            self.system.particles[-1].initial_orientation = [row["initial_ux"], row["initial_uy"], row["initial_uz"]]
            self.system.particles[-1].translation_bound_sq = row["translation_bound_sq"]
            self.system.particles[-1].rotation_bound = row["rotation_bound"]
            self.system.particles[-1].surface_affinity_translation = row["surface_affinity_translation"]
            self.system.particles[-1].surface_affinity_rotation = row["surface_affinity_rotation"]
            self.system.particles[-1].self_affinity = row["self_affinity"]
            self.system.particles[-1].other_affinity = row["other_affinity"]
        
        self.setupCells()
        self.placeParticlesInCells()



    def prepareSimulation(self):
        if not SignalHandler.ProgramState == ProgramState.SHUTDOWN:
            SignalHandler.ProgramState = ProgramState.SETUP

            try:
                input_data_file_path = os.path.join(self.input["dir"], self.input["filename"])
                self.prepareFromData(input_data_file_path)
            except TypeError as e:
                print("TypeError: unable to read input file:", e.__dict__)
                print("creating System new")
                self.prepareNew()
            except OSError as e:
                print("OSError: unable to read datafile:", e.__dict__)
                print("creating System new")
                self.prepareNew()
            # self.makeInteractionLookupTable()
            self.system.makeNeighborLists()



    def sample(self, steps=None, timestats=False, analysis=False):
        # print("sampling")
        if not SignalHandler.ProgramState == ProgramState.SHUTDOWN:
            SignalHandler.ProgramState = ProgramState.RUNNING

            self.sample_starttime = time.perf_counter()

            if not isinstance(steps, type(None)):
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
                self.writeTrajectoryHDF(timestats=timestats, analysis=analysis)
                print(f"epot = {self.system.potentialEnergyConcurrent():.2f}", end="")
                if timestats: print()
                if SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                    print("shutting down")
                    sys.exit(SignalHandler.signal)

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
                        self.writeTrajectoryHDF(timestats=timestats, analysis=analysis)
                        print(f"epot = {self.system.potentialEnergyConcurrent():.2f}", end="")
                        if timestats: print()

                    if SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                        print("shutting down")
                        sys.exit(SignalHandler.signal)
                    self.sample_starttime = time.perf_counter()



    def setFromMetadata(self, metadata):
        self.system.box.x = metadata["box.x"]
        self.system.box.y = metadata["box.y"]
        self.system.box.z = metadata["box.z"]
        self.system.temperature = metadata["temperature"]
        self.system.exchange_global_number = metadata["exchange.global.number"]
        self.system.exchange_global_etot_theshold = metadata["exchange.global.etot_theshold"]
        self.system.exchange_global_orientation = metadata["exchange.global.orientation"]
        self.system.exchange_local_number = metadata["exchange.local.number"]
        self.system.exchange_local_etot_theshold = metadata["exchange.local.etot_theshold"]
        self.system.exchange_local_orientation = metadata["exchange.local.orientation"]
        self.time_actual = metadata["time"]



    def getMetadata(self):
        metadata = {}
        metadata["box.x"] = self.system.box.x
        metadata["box.y"] = self.system.box.y
        metadata["box.z"] = self.system.box.z
        metadata["threads"] = self.system.threads
        metadata["temperature"] = self.system.temperature
        metadata["exchange.global.number"] = self.system.exchange_global_number
        metadata["exchange.global.etot_theshold"] = self.system.exchange_global_etot_theshold
        metadata["exchange.global.orientation"] = self.system.exchange_global_orientation
        metadata["exchange.local.number"] = self.system.exchange_local_number
        metadata["exchange.local.etot_theshold"] = self.system.exchange_local_etot_theshold
        metadata["exchange.local.orientation"] = self.system.exchange_local_orientation
        metadata["translation_step"] = self.system.translation()
        metadata["rotation_step"] = self.system.rotation()
        metadata["interaction.cutoff"] = self.system.interaction_cutoff
        metadata["interaction.surface"] = int(self.system.interaction_surface)
        metadata["interaction.surface_width"] = self.system.interaction_surface_width
        metadata["time"] = self.time_actual
        metadata["cell_min_size"] = self.prms_complete["control"]["cell_min_size"]
        try:
            metadata["simulation_walltime"] = time.perf_counter() - self.sample_starttime
        except:
            metadata["simulation_walltime"] = time.perf_counter()
        return metadata



    def writeTrajectoryHDF(self, log=False, timestats=False, analysis=False):
        if self.output["filename"] == None:
            return

        starttime = time.perf_counter()

        positions = np.array([self.system.box.scaleToBox(p.position) for p in self.system.particles])
        initial_positions = np.array([self.system.box.scaleToBox(p.initial_position) for p in self.system.particles])
        initial_orientations = np.array([self.system.box.scaleToBox(p.initial_orientation) for p in self.system.particles])

        df = pd.DataFrame(dict(
            name = [p.name for p in self.system.particles],
            x = np.array([p[0] for p in positions], dtype=np.float32),
            y = np.array([p[1] for p in positions], dtype=np.float32),
            z = np.array([p[2] for p in positions], dtype=np.float32),
            ux = np.array([p.ux for p in self.system.particles], dtype=np.float32),
            uy = np.array([p.uy for p in self.system.particles], dtype=np.float32),
            uz = np.array([p.uz for p in self.system.particles], dtype=np.float32),
            initial_x = np.array([p[0] for p in initial_positions], dtype=np.float32),
            initial_y = np.array([p[1] for p in initial_positions], dtype=np.float32),
            initial_z = np.array([p[2] for p in initial_positions], dtype=np.float32),
            translation_bound_sq = np.array([p.translation_bound_sq for p in self.system.particles], dtype=np.float32),
            initial_ux = np.array([p[0] for p in initial_orientations], dtype=np.float32),
            initial_uy = np.array([p[1] for p in initial_orientations], dtype=np.float32),
            initial_uz = np.array([p[2] for p in initial_orientations], dtype=np.float32),
            rotation_bound = np.array([p.rotation_bound for p in self.system.particles], dtype=np.float32),
            sigma = np.array([p.sigma for p in self.system.particles], dtype=np.float32),
            epsilon = np.array([p.epsilon for p in self.system.particles], dtype=np.float32),
            kappa = np.array([p.kappa for p in self.system.particles], dtype=np.float32),
            gamma = np.array([p.gamma for p in self.system.particles], dtype=np.float32),
            surface_affinity_translation = np.array([p.surface_affinity_translation for p in self.system.particles], dtype=np.float32),
            surface_affinity_rotation = np.array([p.surface_affinity_rotation for p in self.system.particles], dtype=np.float32),
            self_affinity = np.array([p.self_affinity for p in self.system.particles], dtype=np.float32),
            other_affinity = np.array([p.other_affinity for p in self.system.particles], dtype=np.float32)
        ))

        if analysis or self.direct_analysis:
            if timestats:
                analysis_starttime = time.perf_counter()
            df = analyzeSnapshot(
                df = df, 
                prms = self.prms_complete, 
                metadata = self.getMetadata(),
                system = self.system
            )
            if timestats:
                analysis_endtime = time.perf_counter()

        h5store(
            os.path.join(self.output["dir"], self.output["filename"]),
            f"/time{self.time_actual}",
            df,
            **self.getMetadata()
        )
        
        # df.to_hdf(
        #     path_or_buf=os.path.join(self.output["dir"], self.output["filename"]),
        #     key=f"/time{self.time_actual}",
        #     mode="a",# if self.output["mode"] == "append" else "w",
        #     # append=True if self.output["mode"] == "append" else False,
        #     format="table",
        #     complevel=1,
        # )


        # h5pyFile = h5py.File(os.path.join(self.output["dir"], self.output["filename"]), mode="a")
        # node = h5pyFile.get(f"time{self.time_actual}")
        # node.attrs["box.x"] = self.system.box.x
        # node.attrs["box.y"] = self.system.box.y
        # node.attrs["box.z"] = self.system.box.z
        # node.attrs["threads"] = self.system.threads
        # node.attrs["temperature"] = self.system.temperature
        # node.attrs["translation_step"] = self.system.translation()
        # node.attrs["rotation_step"] = self.system.rotation()
        # node.attrs["interaction.cutoff"] = self.system.interaction_cutoff
        # node.attrs["interaction_surface"] = int(self.system.interaction_surface)
        # node.attrs["interaction_surface_width"] = self.system.interaction_surface_width
        # node.attrs["time"] = self.time_actual
        # node.attrs["cell_min_size"] = self.prms_complete["control"]["cell_min_size"]
        # try:
        #     node.attrs["simulation_walltime"] = time.perf_counter() - self.sample_starttime
        # except:
        #     node.attrs["simulation_walltime"] = time.perf_counter()

        # h5pyFile.close()

        if log:
            print(df)
            print(df.info())
            
        if timestats:
            endtime = time.perf_counter()
            if analysis or self.direct_analysis:
                print(f" write HDF5 took {endtime-starttime:.4f} s  ({(analysis_endtime-analysis_starttime)/(endtime-starttime)*100:.0f}% analysis) |  ", end="")
            else:
                print(f" write HDF5 took {endtime-starttime:.4f} s  |  ", end="")