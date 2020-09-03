import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_construction(self):
        sys = pyves.System()
        print(sys.particles)



if __name__ == '__main__':
    unittest.main()