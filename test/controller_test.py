from binascii import a2b_qp
import enum
import unittest
import pyves
import os
import numpy as np



class MainTest(unittest.TestCase):
    def test_simulation_run(self):
        print()
        try:
            os.remove("test/data.h5")
        except Exception as e:
            print("no data.h5 to remove",e)

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

        control.sample(steps=500, timestats=True)

        for cell in control.system.cells:
            for particle in cell.particles:
                self.assertTrue(cell.contains(particle))
        

        control2 = pyves.Controller()
        control2.readForceField("test/forcefield.json")
        control2.readParameters("test/parameters.json")

        control2.prepareSimulation()
        box_dims = np.array([control.system.box.x, control.system.box.y, control.system.box.z])
        cells_per_dim = np.array(box_dims/control.cell_min_size).astype(int)
        self.assertEqual(len(control2.system.cells), np.cumprod(cells_per_dim)[-1])

        control2.sample(timestats=True)
        self.assertEqual(control2.system.numParticlesInCells(), 100)
        for i,p in enumerate(control2.system.particles):
            self.assertTrue(p.assertIntegrity(), f"particle {i}")
        for i,c in enumerate(control2.system.cells):
            for j,p in enumerate(c.particles):
                print()
                print(c, p)
                print(c.contains(p))
                print(c.insideCellBounds(p))
            self.assertTrue(c.assertIntegrity(), f"cell {i}")
        
        self.assertTrue(control2.system.assertIntegrity())
        # import pprint
        # pprint.pprint(control.__dict__)
        # print(pyves.interaction(control.system.particles[0], control.system.particles[1], control.system.box, 3))



if __name__ == '__main__':
    unittest.main()