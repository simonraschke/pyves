import unittest
# import _pyves as pyves
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_epot_x_direction(self):
        p1 = pyves.Particle([0,0,0], np.array([-1,1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([1.5*2**(1.0/6),0,0], np.array([1,1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 2
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)
        
        p1 = pyves.Particle([0,0,0], np.array([-1,-1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([1.5*2**(1.0/6),0,0], np.array([1,-1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 2
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)
        
        p1 = pyves.Particle([0,0,0], [0,1,0])
        p2 = pyves.Particle([2**(1.0/6),0,0], [0,1,0])
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)

        p1 = pyves.Particle([0,0,0], [0,-1,0])
        p2 = pyves.Particle([2**(1.0/6),0,0], [0,-1,0])
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)
        e = pyves.interaction(p1, p2, box,1)
        self.assertAlmostEqual(e, 0)



    def test_epot_y_direction(self):
        p1 = pyves.Particle([0,0,0], np.array([1,-1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([0,1.5*2**(1.0/6),0], np.array([1,1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 2
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)
        
        p1 = pyves.Particle([0,0,0], np.array([-1,-1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([0,1.5*2**(1.0/6),0], np.array([-1,1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 2
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)
        
        p1 = pyves.Particle([0,0,0], [1,0,0])
        p2 = pyves.Particle([0,2**(1.0/6),0], [1,0,0])
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)
        
        p1 = pyves.Particle([0,0,0], [1,0,0])
        p2 = pyves.Particle([0,20 - 2**(1.0/6),0], [1,0,0])
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)

        # p1 = pyves.Particle([0,0,0], [0,-1,0])
        # p2 = pyves.Particle([2**(1.0/6),0,0], [0,-1,0])
        # p1.gamma = np.pi/180*0
        # p2.gamma = np.pi/180*0
        # p1.kappa = 1 
        # p2.kappa = 1
        # p1.epsilon = p2.epsilon = 1
        # p1.sigma = 1
        # p2.sigma = 1
        # box = pyves.BoxPBC(10,10,10)
        # e = pyves.interaction(p1, p2, box,3)
        # self.assertAlmostEqual(e, -1)



    def test_epot_diagonal(self):
        p1 = pyves.Particle([0,0,0], [-1,0,0])
        p2 = pyves.Particle([1,1,0]/np.sqrt(2)*2**(1.0/6), [0,1,0])
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        self.assertAlmostEqual(e, -1)



    def test_epot_translation(self):
        p1 = pyves.Particle([0,0,0], np.array([-1,1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([1,1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        # print( )
        # print(2**(1.0/6))
        # for x in np.arange(1.0, 3.1, 2**(1.0/6)-1):
        #     p2.forceSetPosition([x,0,0])
        #     print(f"x={x:.2f}  epot={pyves.interaction(p1, p2, box,3):.4f}")
        # self.assertAlmostEqual(e, -1)



    def test_epot_rotation(self):
        from scipy.spatial.transform import Rotation as R
        p1 = pyves.Particle([0,0,0], [0,1,0])
        p2 = pyves.Particle([2**(1.0/6),0,0], [0,1,0])
        p1.gamma = np.pi/180*10
        p2.gamma = np.pi/180*10
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box,3)
        # print()
        p1.forceSetOrientation(R.from_euler('z', 10, degrees=True).apply(p1.orientation))
        p2.forceSetOrientation(R.from_euler('z', -10, degrees=True).apply(p2.orientation))
        # print(f"theta={10:.2f}  epot={pyves.interaction(p1, p2, box,3):.4f}")
        # for t in np.arange(10,360,40):
        #     p2.forceSetOrientation(R.from_euler('z', -40, degrees=True).apply(p2.orientation))
        #     print(f"theta={t:.2f}  epot={pyves.interaction(p1, p2, box,3):.4f}")


    
    def test_interaction_integrity(self):
        ctrl = pyves.Controller()
        ctrl.readParameters("test/benchmark.json")
        ctrl.prepareSimulation()
        ctrl.sample()

        r1 = ctrl.system.potentialEnergy()
        r2 = ctrl.system.potentialEnergyConcurrent()

        res = [r1,r2]

        # print("\n\nINTERACTION RESULTS")
        # print( "1 thread   ", r1)
        # print(f"{ctrl.system.threads} threads  ", r2)
        
        for i, resi in enumerate(res):
            for j, resj in enumerate(res):
                if i != j:
                    self.assertAlmostEqual(resi, resj, 1)




    def test_interaction_benchmark(self):
        ctrl = pyves.Controller()
        ctrl.readParameters("test/benchmark.json")
        ctrl.prepareSimulation()
        ctrl.sample()
        
        print("\n\nINTERACTION TIMES")
        ctrl.system.benchmark(100)



    def test_surface_potential(self):

        def angle(v1, v2):
            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)
            return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2)))


        surface_width = 2
        interaction_cutoff = 3
        UNIT_Z = [0,0,1]

        p = pyves.Particle([1,1,19.999999], [0,0,1])
        p.surface_affinity_translation = 1
        p.surface_affinity_rotation = 1

        box = pyves.BoxPBC(10,10,20)

        _a = pyves.internal_angle_pow2_penalty(p)
        self.assertAlmostEqual(_a, 0, 5)
        
        o = [0,0,-1]
        p.forceSetOrientation(o)
        _a = pyves.internal_angle_pow2_penalty(p)
        self.assertAlmostEqual(_a, np.pi*np.pi, 5)
        self.assertAlmostEqual(_a, angle(p.orientation, UNIT_Z)*angle(p.orientation, UNIT_Z), 5)
        
        o = [1,0,1]
        p.forceSetOrientation(o)
        _a = pyves.internal_angle_pow2_penalty(p)
        self.assertTrue(_a > 0)
        self.assertAlmostEqual(_a, angle(p.orientation, UNIT_Z)*angle(p.orientation, UNIT_Z), 5)

        o = [-1,-1,1]
        p.forceSetOrientation(o)
        _a = pyves.internal_angle_pow2_penalty(p)
        self.assertTrue(_a > 0)
        self.assertAlmostEqual(_a, angle(p.orientation, UNIT_Z)*angle(p.orientation, UNIT_Z), 5)

        o = [1,0,0]
        p.forceSetOrientation(o)
        _a = pyves.internal_angle_pow2_penalty(p)
        self.assertAlmostEqual(_a, np.pi*np.pi/4, 5)
        self.assertAlmostEqual(_a, angle(p.orientation, UNIT_Z)*angle(p.orientation, UNIT_Z), 5)

        p.forceSetOrientation(UNIT_Z)
        value = pyves.internal_z_direction_energy(p.z, box, surface_width, interaction_cutoff)
        self.assertAlmostEqual(value, -(p.z- (box.z-surface_width))/surface_width, 5)
        self.assertAlmostEqual(value, -1, 5)

        surface_width = 2
        p.forceSetPosition([1,1,19])
        p.forceSetOrientation(UNIT_Z)
        value = pyves.internal_z_direction_energy(p.z, box, surface_width, interaction_cutoff)
        self.assertAlmostEqual(value, -(p.z- (box.z-surface_width))/surface_width, 5)
        self.assertAlmostEqual(value, -0.5, 5)

        surface_width = 1
        p.forceSetOrientation(UNIT_Z)
        value = pyves.internal_z_direction_energy(p.z, box, surface_width, interaction_cutoff)
        self.assertAlmostEqual(value, -(p.z- (box.z-surface_width))/surface_width, 5)
        self.assertAlmostEqual(value, 0, 5)

        p.forceSetPosition([1,1,19.5])
        p.forceSetOrientation(UNIT_Z)
        value = pyves.internal_z_direction_energy(p.z, box, surface_width, interaction_cutoff)
        self.assertAlmostEqual(value, -(p.z- (box.z-surface_width))/surface_width, 5)
        self.assertAlmostEqual(value, -0.5, 5)

        p.forceSetPosition([1,1,19.75])
        p.forceSetOrientation(UNIT_Z)
        value = pyves.internal_z_direction_energy(p.z, box, surface_width, interaction_cutoff)
        self.assertAlmostEqual(value, -(p.z- (box.z-surface_width))/surface_width, 5)
        self.assertAlmostEqual(value, -0.75, 5)
        self.assertAlmostEqual(value, pyves.surface_potential(p, box, surface_width, interaction_cutoff), 5)
        self.assertAlmostEqual(value, pyves.external_potential(p, box, surface_width, interaction_cutoff), 5)
        
        p.surface_affinity_translation = 0.2
        p.forceSetPosition([1,1,19.999999])
        self.assertAlmostEqual(-0.2, pyves.surface_potential(p, box, surface_width, interaction_cutoff), 5)
        self.assertAlmostEqual(-0.2, pyves.external_potential(p, box, surface_width, interaction_cutoff), 5)
        
        surface_width = 2
        p.surface_affinity_translation = 1
        p.forceSetPosition([1,1,2.99])
        self.assertTrue(pyves.surface_potential(p, box, surface_width, interaction_cutoff) > 1e5 )
        self.assertTrue(pyves.external_potential(p, box, surface_width, interaction_cutoff) > 1e5 )

        p.forceSetPosition([1,1,18.5])
        self.assertTrue(pyves.surface_potential(p, box, surface_width, interaction_cutoff) < 1e5 )
        self.assertTrue(pyves.external_potential(p, box, surface_width, interaction_cutoff) < 1e5 )

        p.forceSetPosition([1,1,interaction_cutoff+0.001])
        self.assertAlmostEqual(pyves.surface_potential(p, box, surface_width, interaction_cutoff), 0)
        self.assertAlmostEqual(pyves.external_potential(p, box, surface_width, interaction_cutoff), 0 )

        p.forceSetPosition([1,1,interaction_cutoff-0.001])
        self.assertTrue(pyves.surface_potential(p, box, surface_width, interaction_cutoff) > 1e5 )
        self.assertTrue(pyves.external_potential(p, box, surface_width, interaction_cutoff) > 1e5 )

        p.forceSetPosition([1,1,.5])
        self.assertTrue(pyves.surface_potential(p, box, surface_width, interaction_cutoff) > 1e5 )
        self.assertTrue(pyves.external_potential(p, box, surface_width, interaction_cutoff) > 1e5 )



    def test_other_affinity(self):
        p1 = pyves.Particle([0,0,0], np.array([-1,1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([1,1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = 1
        p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        p1.name = "A"
        p2.name = "B"
        p1.self_affinity = 1
        p2.self_affinity = 1
        p1.other_affinity = 0.9
        p2.other_affinity = 0.9
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box, 3)
        self.assertAlmostEqual(e, -0.9)

        p1 = pyves.Particle([0,0,0], np.array([-1,1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([1,1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = 1
        p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        p1.name = "A"
        p2.name = "A"
        p1.self_affinity = 1
        p2.self_affinity = 1
        p1.other_affinity = 0.9
        p2.other_affinity = 0.9
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box, 3)
        self.assertAlmostEqual(e, -1)

        p1 = pyves.Particle([0,0,0], np.array([0,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([0,1,0]))
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = 1
        p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        p1.name = "A"
        p2.name = "B"
        p1.self_affinity = 1
        p2.self_affinity = 1
        p1.other_affinity = 0.7
        p2.other_affinity = 0.9
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box, 3)
        self.assertAlmostEqual(e, -0.8)

        p1 = pyves.Particle([0,0,0], np.array([0,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([0,1,0]))
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = 1
        p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        p1.name = "A"
        p2.name = "B"
        p1.self_affinity = 1
        p2.self_affinity = 1
        p1.other_affinity = 1.7
        p2.other_affinity = 1.9
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box, 3)
        self.assertAlmostEqual(e, -1.8)



    def test_self_affinity(self):
        p1 = pyves.Particle([0,0,0], np.array([-1,1,0])/np.linalg.norm([-1,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([1,1,0])/np.linalg.norm([1,1,0]))
        p1.gamma = np.pi/180*45
        p2.gamma = np.pi/180*45
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = 1
        p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        p1.name = "A"
        p2.name = "A"
        p1.self_affinity = 1.1
        p2.self_affinity = 1.1
        p1.other_affinity = 1
        p2.other_affinity = 1
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box, 3)
        self.assertAlmostEqual(e, -1.1)

        p1 = pyves.Particle([0,0,0], np.array([0,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([0,1,0]))
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = 1
        p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        p1.name = "A"
        p2.name = "A"
        p1.self_affinity = 1
        p2.self_affinity = 1
        p1.other_affinity = 0.7
        p2.other_affinity = 0.9
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box, 3)
        self.assertAlmostEqual(e, -1)

        p1 = pyves.Particle([0,0,0], np.array([0,1,0]))
        p2 = pyves.Particle([2**(1.0/6),0,0], np.array([0,1,0]))
        p1.gamma = np.pi/180*0
        p2.gamma = np.pi/180*0
        p1.kappa = 1 
        p2.kappa = 1
        p1.epsilon = 1
        p2.epsilon = 1
        p1.sigma = 1
        p2.sigma = 1
        p1.name = "A"
        p2.name = "A"
        p1.self_affinity = 1.8
        p2.self_affinity = 1.8
        p1.other_affinity = 0.7
        p2.other_affinity = 0.9
        box = pyves.BoxPBC(10,10,10)
        e = pyves.interaction(p1, p2, box, 3)
        self.assertAlmostEqual(e, -1.8)



if __name__ == '__main__':
    unittest.main()
