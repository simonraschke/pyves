import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_construction(self):
        p1 = pyves.Particle()
        p2 = pyves.Particle(p1)
        p3 = pyves.Particle([1,1,1], np.array([1,0,0]))

    def test_operators(self):
        p1 = pyves.Particle()
        p2 = pyves.Particle()
        p3 = pyves.Particle(p1)
        self.assertTrue(p1 == p1)
        self.assertTrue(p2 == p2)
        self.assertFalse(p2 == p3)
        self.assertFalse(p1 == p3)
        self.assertFalse(p1 != p1)
        self.assertFalse(p2 != p2)
        self.assertTrue(p2 != p3)
        self.assertTrue(p1 != p3)

        p4 = pyves.Particle([1,1,1], np.array([1,0,0]))
        p5 = pyves.Particle()
        p5.position = p4.position + p3.position
        self.assertAlmostEqual(p5.position[0], 1)
        self.assertAlmostEqual(p5.position[1], 1)
        self.assertAlmostEqual(p5.position[2], 1)
        self.assertAlmostEqual(p5.x(), 1)
        self.assertAlmostEqual(p5.y(), 1)
        self.assertAlmostEqual(p5.z(), 1)
        self.assertAlmostEqual(p5.position[0], p5.x())
        self.assertAlmostEqual(p5.position[1], p5.y())
        self.assertAlmostEqual(p5.position[2], p5.z())


    # def test_operator_plus(self):





if __name__ == '__main__':
    unittest.main()