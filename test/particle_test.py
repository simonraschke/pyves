import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_construction(self):
        p1 = pyves.Particle()
        p1.position = [1,p1.y,p1.z]
        self.assertAlmostEqual(p1.x, 1)

        p2 = pyves.Particle(p1)
        self.assertAlmostEqual(p2.x, 1)

        p2.position = np.array([2,p2.y,p2.z])
        self.assertAlmostEqual(p1.x, 1)
        self.assertAlmostEqual(p2.x, 2)

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
        self.assertAlmostEqual(p5.x, 1)
        self.assertAlmostEqual(p5.y, 1)
        self.assertAlmostEqual(p5.z, 1)
        self.assertAlmostEqual(p5.position[0], p5.x)
        self.assertAlmostEqual(p5.position[1], p5.y)
        self.assertAlmostEqual(p5.position[2], p5.z)
        self.assertAlmostEqual(p3.position[0], 0)
        self.assertAlmostEqual(p3.position[1], 0)
        self.assertAlmostEqual(p3.position[2], 0)
        self.assertAlmostEqual(p3.x, 0)
        self.assertAlmostEqual(p3.y, 0)
        self.assertAlmostEqual(p3.z, 0)

        p6 = pyves.Particle()
        p6.position = [1,2,3]
        p6.orientation = [45.1239874,0,0]
        self.assertAlmostEqual(p6.position[0], 1)
        self.assertAlmostEqual(p6.position[1], 2)
        self.assertAlmostEqual(p6.position[2], 3)
        self.assertAlmostEqual(p6.orientation[0], 1)
        self.assertAlmostEqual(p6.orientation[1], 0)
        self.assertAlmostEqual(p6.orientation[2], 0)


    # def test_operator_plus(self):





if __name__ == '__main__':
    unittest.main()