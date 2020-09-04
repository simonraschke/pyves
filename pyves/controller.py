import json
import _pyves
import numpy as np



class Controller(object):
    def __init__(self) -> None:
        self.system = _pyves.System()
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
            self.system.translation_min = _prms["system"]["translation_min"]
            self.system.translation_max = _prms["system"]["translation_max"]
            self.system.rotation_min = _prms["system"]["rotation_min"]
            self.system.rotation_max = _prms["system"]["rotation_max"]
            self.system.time_max = _prms["system"]["time_max"]

            self.time_delta = _prms["control"]["time_delta"]
            self.num_particles = _prms["control"]["num_particles"]
            self.cell_min_size = _prms["control"]["cell_min_size"]



    def setParameters(self, prms:dict):
        self.parameters = prms
    


    def prepareSimulation(self):
        # Setup particles
        assert( hasattr(self, "forcefield") )
        for name, num in self.num_particles.items():
            assert( name in self.forcefield )
            pprms = self.forcefield[name]
            for _ in range(num):
                self.system.particles.append(_pyves.Particle([0,0,0], [1,0,0], 
                    sigma=pprms["sigma"], kappa=pprms["kappa"], eps=pprms["epsilon"], gamma=pprms["gamma"], name=name))
        
        # Setup cells
        box_dims = np.array([self.system.box.x, self.system.box.y, self.system.box.z])
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
                if ci.isNeighbourOf(cj, self.system.box):
                    ci.proximity.append(cj)
        for cell in self.system.cells:
            assert( cell.assertIntegrity() )
        
        assert(self.system.assertIntegrity())