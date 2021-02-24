import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    
    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(MainTest, self).run(result)


    def test_construction(self):
        sys = pyves.System()



    def test_particle_exchange(self):
        sys = pyves.System()
        sys.box = pyves.BoxPBC(10,10,10)
        sys.particles.append(pyves.Particle([1,1,1], [1,0,0], sigma=2, kappa=1, eps=4, name="A", gamma=0))        
        sys.particles.append(pyves.Particle([3,1,1], [1,0,0], sigma=1, kappa=3, eps=1, name="B", gamma=0))

        self.assertEqual(sys.particles[0].name, "A")
        self.assertEqual(sys.particles[1].name, "B")
        print()
        print(sys.particles[0].detailed_repr())
        print(sys.particles[1].detailed_repr())

        sys.exchangeParticleParameters(sys.particles[0], sys.particles[1])
        print(sys.particles[0].detailed_repr())
        print(sys.particles[1].detailed_repr())

        self.assertEqual(sys.particles[0].name, "B")
        self.assertEqual(sys.particles[1].name, "A")
        self.assertAlmostEqual(sys.particles[0].sigma, 1, 1e6)
        self.assertAlmostEqual(sys.particles[1].sigma, 2, 1e6)
        self.assertAlmostEqual(sys.particles[0].kappa, 3, 1e6)
        self.assertAlmostEqual(sys.particles[1].kappa, 1, 1e6)
        # self.assertEqual(sys.particles[0].position[0], 3)
        # self.assertEqual(sys.particles[1].position[0], 1)



    def test_global_exchange(self):
        sys = pyves.System()
        sys.temperature = 0.3
        sys.exchange_global_number = 50
        sys.exchange_global_etot_theshold = 100000
        sys.box = pyves.BoxPBC(10,10,10)
        sys.particles.append(pyves.Particle([1,1,1], [1,1,0], sigma=1, kappa=1, eps=1, name="A", gamma=np.pi/5))        
        sys.particles.append(pyves.Particle([2.1224,1,1], [1,-1,0], sigma=1, kappa=1, eps=1, name="B", gamma=np.pi/4))

        sys.prepareSimulationStep()

        self.assertEqual(sys.particles[0].name, "A")
        self.assertEqual(sys.particles[1].name, "B")
        self.assertAlmostEqual(sys.particles[0].position[0], 1.0, 5)
        self.assertAlmostEqual(sys.particles[0].position[1], 1.0, 5)
        self.assertAlmostEqual(sys.particles[0].position[2], 1.0, 5)
        self.assertAlmostEqual(sys.particles[0].orientation[0], 1.0/np.linalg.norm([1,1,0]), 5)
        self.assertAlmostEqual(sys.particles[0].orientation[1], 1.0/np.linalg.norm([1,1,0]), 5)
        self.assertAlmostEqual(sys.particles[0].orientation[2], 0.0/np.linalg.norm([1,1,0]), 5)

        print("\n------------------------")
        print("sys.exchange_global_number", sys.exchange_global_number)
        print("sys.particles.size", len(sys.particles))
        print("------------------------")
        sys.exchangeParticleParameters(sys.particles[0], sys.particles[1])

        self.assertAlmostEqual(sys.particles[1].x, 2.1224, 5)
        self.assertAlmostEqual(sys.particles[1].y, 1.0, 5)
        self.assertAlmostEqual(sys.particles[1].z, 1.0, 5)
        self.assertAlmostEqual(sys.particles[1].ux, 1.0/np.linalg.norm([1,1,0]), 5)
        self.assertAlmostEqual(sys.particles[1].uy, -1.0/np.linalg.norm([1,1,0]), 5)
        self.assertAlmostEqual(sys.particles[1].uz, 0.0/np.linalg.norm([1,1,0]), 5)
        self.assertEqual(sys.particles[0].name, "B")
        self.assertEqual(sys.particles[1].name, "A")
        self.assertAlmostEqual(sys.particles[0].sigma, 1.0, 5)
        self.assertAlmostEqual(sys.particles[1].sigma, 1.0, 5)
        self.assertAlmostEqual(sys.particles[0].kappa, 1.0, 5)
        self.assertAlmostEqual(sys.particles[1].kappa, 1.0, 5)
        self.assertAlmostEqual(sys.particles[0].epsilon, 1.0, 5)
        self.assertAlmostEqual(sys.particles[1].epsilon, 1.0, 5)
        self.assertAlmostEqual(sys.particles[0].gamma, np.pi/4, 5)
        self.assertAlmostEqual(sys.particles[1].gamma, np.pi/5, 5)
        
        # np.array(a)
        # self.assertEqual(sys.particles[0].name, "B")
        # self.assertEqual(sys.particles[1].name, "A")

        # ctrl = pyves.Controller()
        # ctrl.readParameters("test/parameters.json")
        # ctrl.prepareSimulation()
        # ctrl.system.particles = [
        #     pyves.Particle([1,1,1], [1,.5,0], sigma=1, kappa=1, eps=1, name="A", gamma=0),
        # ]
        # ctrl.system.temperature = 0.3
        # ctrl.system.exchange_global_ratio = 1.0
        # ctrl.system.box = pyves.BoxPBC(10,10,10)
        # ctrl.system.particles.append(pyves.Particle([1,1,1], [1,.5,0], sigma=1, kappa=1, eps=1, name="A", gamma=0))
        # ctrl.system.particles.append(pyves.Particle([2.1224,1,1], [1,-0.5,0], sigma=1, kappa=1, eps=1, name="B", gamma=0))

        # self.assertEqual(ctrl.system.particles[0].name, "A")
        # self.assertEqual(ctrl.system.particles[1].name, "B")

        # ctrl.system.globalExchange(1e5)

        # self.assertEqual(ctrl.system.particles[0].name, "B")
        # self.assertEqual(ctrl.system.particles[1].name, "A")



if __name__ == '__main__':
    unittest.main(failfast=True)