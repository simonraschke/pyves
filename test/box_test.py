import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_construction(self):
        box1 = pyves.BoxPBC()
        self.assertAlmostEqual(box1.x, 0)
        self.assertAlmostEqual(box1.y, 0)
        self.assertAlmostEqual(box1.z, 0)

        box2 = pyves.BoxPBC(20,10,15)
        self.assertAlmostEqual(box2.x, 20)
        self.assertAlmostEqual(box2.y, 10)
        self.assertAlmostEqual(box2.z, 15)

        box3 = pyves.BoxPBC([10,15,20])
        self.assertAlmostEqual(box3.x, 10)
        self.assertAlmostEqual(box3.y, 15)
        self.assertAlmostEqual(box3.z, 20)

        box4 = pyves.BoxPBC(np.array([10,15,20]))
        self.assertAlmostEqual(box4.x, 10)
        self.assertAlmostEqual(box4.y, 15)
        self.assertAlmostEqual(box4.z, 20)



    def test_simple_methods(self):
        box1 = pyves.BoxPBC()
        box2 = pyves.BoxPBC(20,10,15)
        box3 = pyves.BoxPBC([10,15,20])
        box4 = pyves.BoxPBC(np.array([10,15,20]))

        # print()
        # print(pyves.BoxPBC.setx.__dict__)
        box1.x = 7
        self.assertAlmostEqual(box1.x, 7)

        box1.y = 3
        self.assertAlmostEqual(box1.y, 3)

        box1.z = 6
        self.assertAlmostEqual(box1.z, 6)

        self.assertAlmostEqual(box1.volume, 7*3*6)
        self.assertAlmostEqual(box2.volume, 20*10*15)
        self.assertAlmostEqual(box3.volume, 10*15*20)

        self.assertAlmostEqual(np.linalg.norm(box2.center-np.array([10,5,7.5])), 0)



    def test_PBC_NoPBC_methods(self):
        boxPBC = pyves.BoxPBC(10,10,10)
        boxNoPBC = pyves.BoxNoPBC([10,10,10])

        a000 = np.array([0,0,0])
        a111 = np.array([1,1,1])
        b111 = np.array([11,11,11])

        self.assertTrue(np.allclose(boxPBC.distanceVector(a000, a111), a111.copy()))
        self.assertTrue(np.allclose(boxPBC.distanceVector(a000, b111), a111.copy()))
        self.assertTrue(np.allclose(boxNoPBC.distanceVector(a000, a111), a111.copy()))
        self.assertTrue(np.allclose(boxNoPBC.distanceVector(a000, b111), b111.copy()))

        self.assertAlmostEqual(boxPBC.distance(a000, a111), 1.7320508075688772)
        self.assertAlmostEqual(boxPBC.distance(a000, b111), 1.7320508075688772)
        self.assertAlmostEqual(boxNoPBC.distance(a000, a111), 1.7320508075688772)
        self.assertAlmostEqual(boxNoPBC.distance(a000, b111), 19.05255888325765)

        self.assertAlmostEqual(boxPBC.squaredDistance(a000, a111), 3)
        self.assertAlmostEqual(boxPBC.squaredDistance(a000, b111), 3)
        self.assertAlmostEqual(boxNoPBC.squaredDistance(a000, a111), 3)
        self.assertAlmostEqual(boxNoPBC.squaredDistance(a000, b111), 363)

        self.assertTrue(np.allclose(boxPBC.scaleToBox(b111.copy()), a111.copy()))

        self.assertTrue(boxPBC.contains(np.array([4,4,4])))
        self.assertTrue(boxPBC.contains(np.array([10,10,10])))
        self.assertTrue(boxPBC.contains(np.array([11,11,11])))
        self.assertTrue(boxPBC.contains(np.array([1,-1,-21])))
        self.assertTrue(boxNoPBC.contains(np.array([4,4,4])))
        self.assertTrue(boxNoPBC.contains(np.array([10,10,10])))
        self.assertFalse(boxNoPBC.contains(np.array([11,11,11])))
        self.assertFalse(boxNoPBC.contains(np.array([1,-1,-21])))

        self.assertTrue(boxNoPBC.contains(boxNoPBC.randomPointInside()))
        self.assertIsInstance(boxNoPBC.randomPointInside(), np.ndarray)
        


if __name__ == '__main__':
    unittest.main()