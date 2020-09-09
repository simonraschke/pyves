import unittest

from pandas.core.resample import f
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_simulation_run(self):
        control = pyves.Controller()
        control.readForceField("test/forcefield.json")
        self.assertAlmostEqual(control.forcefield["A"]["sigma"], 1)

        copy = control.forcefield.copy()
        copy["A"]["sigma"] = 2
        self.assertAlmostEqual(copy["A"]["sigma"], 2)

        control.setForceField(copy)
        self.assertAlmostEqual(control.forcefield["A"]["sigma"], 2)
        control.forcefield["A"]["sigma"] = 1
        
        control.readParameters("test/parameters.json")

        # print(control.forcefield)
        control.prepareSimulation()
        self.assertEqual(len(control.system.particles), 100)
        box_dims = np.array([control.system.box.x, control.system.box.y, control.system.box.z])
        cells_per_dim = np.array(box_dims/control.cell_min_size).astype(int)
        self.assertEqual(len(control.system.cells), np.cumprod(cells_per_dim)[-1])

        self.assertEqual(control.system.numParticlesInCells(), 100)
        self.assertTrue(control.system.assertIntegrity())

        control.sample()

        for cell in control.system.cells:
            for particle in cell.particles:
                self.assertTrue(cell.contains(particle))
        
        # import pprint
        # pprint.pprint(control.__dict__)
        # print(pyves.interaction(control.system.particles[0], control.system.particles[1], control.system.box, 3))



if __name__ == '__main__':
    unittest.main()