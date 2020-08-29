import unittest
import pyves

class MainTest(unittest.TestCase):
    def test_operator_plus(self):
        p1 = pyves.Pointi(1,1,1)
        p2 = pyves.Pointi(4,2,3)
        p3 = p1 + p2
        self.assertEqual(p3.x, 5)
        self.assertEqual(p3.y, 3)
        self.assertEqual(p3.z, 4)

        p1 = pyves.Pointf(1.2,1.3,1.4)
        p2 = pyves.Pointf(4,2,3)
        p3 = p1 + p2
        self.assertAlmostEqual(p3.x, 5.2, 6)
        self.assertAlmostEqual(p3.y, 3.3, 6)
        self.assertAlmostEqual(p3.z, 4.4, 6)

        p1 = pyves.Pointd(1.2,1.3,1.4)
        p2 = pyves.Pointd(4,2,3)
        p3 = p1 + p2
        self.assertAlmostEqual(p3.x, 5.2)
        self.assertAlmostEqual(p3.y, 3.3)
        self.assertAlmostEqual(p3.z, 4.4)



    def test_operator_minus(self):
        p1 = pyves.Pointi(1,1,1)
        p2 = pyves.Pointi(4,2,3)
        p3 = p2-p1
        self.assertEqual(p3.x, 3)
        self.assertEqual(p3.y, 1)
        self.assertEqual(p3.z, 2)

        p1 = pyves.Pointf(1.2,1.3,1.4)
        p2 = pyves.Pointf(4,2,3)
        p3 = p2-p1
        self.assertAlmostEqual(p3.x, 2.8, 6)
        self.assertAlmostEqual(p3.y, 0.7, 6)
        self.assertAlmostEqual(p3.z, 1.6, 6)

        p1 = pyves.Pointd(1.2,1.3,1.4)
        p2 = pyves.Pointd(4,2,3)
        p3 = p2-p1
        self.assertAlmostEqual(p3.x, 2.8)
        self.assertAlmostEqual(p3.y, 0.7)
        self.assertAlmostEqual(p3.z, 1.6)



    def test_oprator_dot(self):
        p1 = pyves.Pointi(1,1,1)
        p2 = pyves.Pointi(4,2,3)
        self.assertEqual(p1.dot(p2), 9)

        p1 = pyves.Pointf(1.2,1.3,1.4)
        p2 = pyves.Pointf(4,2,3)
        self.assertAlmostEqual(p1.dot(p2), 1.2*4 + 1.3*2 + 1.4*3, 6)
        self.assertAlmostEqual(p2.dot(p1), 1.2*4 + 1.3*2 + 1.4*3, 6)

        p1 = pyves.Pointd(1.2,1.3,1.4)
        p2 = pyves.Pointd(4,2,3)
        self.assertAlmostEqual(p1.dot(p2), 1.2*4 + 1.3*2 + 1.4*3)
        self.assertAlmostEqual(p2.dot(p1), 1.2*4 + 1.3*2 + 1.4*3)


if __name__ == '__main__':
    unittest.main()
