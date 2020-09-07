import unittest
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
        self.assertEqual(len(control.system.particles), 200)
        self.assertEqual(len(control.system.cells), 5*6*7)

        control.sample(3)



if __name__ == '__main__':
    unittest.main()