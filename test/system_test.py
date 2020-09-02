import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_construction(self):
        sys = pyves.System()
        print()
        print(sys.particles)
        sys.particles.append(pyves.Particle([1,1,1], [1,0,0]))
        print(sys.particles)
        print()
        print()



if __name__ == '__main__':
    unittest.main()