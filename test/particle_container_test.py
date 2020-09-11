import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_particle_container(self):
        cont = pyves.ParticleContainer()
        self.assertEqual(len(cont), 0)
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        cont.append(pyves.Particle([1,1,1], [1,0,0]))
        self.assertEqual(len(cont), 4)
        cont[0].position = [2, cont[0].y, cont[0].z]
        self.assertAlmostEqual(cont[0].x, 2)



    def test_system_integrity(self):
        sys = pyves.System()
        self.assertEqual(len(sys.particles), 0)
        sys.particles.append(pyves.Particle([1,1,1], [1,0,0], sigma=1, kappa=1, eps=1, name="TEST", gamma=0))
        sys.particles.append(pyves.Particle([1,1,1], [1,0,0], sigma=1, kappa=1, eps=1, name="TEST", gamma=0))
        try:
            sys.particles.append(pyves.Particle([1,1,1], [1,0,0], sigma=np.nan, kappa=1, eps=1, name="TEST", gamma=0))
        except RuntimeError:
            pass
        self.assertEqual(len(sys.particles), 3)
        # sys.particles[0].x = 2
        sys.particles[0].position = [2, sys.particles[0].y, sys.particles[0].z]
        self.assertAlmostEqual(sys.particles[0].x, 2)
        self.assertFalse(sys.assertIntegrity())
        sys.particles[2].sigma = 1
        for particle in sys.particles:
            self.assertTrue(particle.assertIntegrity())



if __name__ == '__main__':
    unittest.main()