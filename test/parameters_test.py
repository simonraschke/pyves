import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_construction(self):
        prms = pyves.Parameters()
        self.assertTrue(np.isnan(prms.system.temperature))

        prms.system.temperature = 0.22
        self.assertAlmostEqual(prms.system.temperature, 0.22)

        prms2 = pyves.Parameters(prms)
        self.assertAlmostEqual(prms2.system.temperature, 0.22)

        prms2.system.temperature = 0.33
        self.assertAlmostEqual(prms.system.temperature, 0.22)
        self.assertAlmostEqual(prms2.system.temperature, 0.33)

        temperature = prms2.system.temperature
        self.assertAlmostEqual(prms2.system.temperature, 0.33)

        temperature = 0.44
        self.assertAlmostEqual(prms.system.temperature, 0.22)
        self.assertAlmostEqual(prms2.system.temperature, 0.33)
        self.assertNotAlmostEqual(prms2.system.temperature, 0.44)



if __name__ == '__main__':
    unittest.main()