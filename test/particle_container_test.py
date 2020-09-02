import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_construction(self):
        cont = pyves.ParticleContainer()
        self.assertEqual(len(cont), 0)
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        self.assertEqual(len(cont), 4)
        cont[0].x = 2
        self.assertAlmostEqual(cont[0].x, 2)



if __name__ == '__main__':
    unittest.main()