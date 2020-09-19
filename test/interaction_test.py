import unittest
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
        
        import pyves
        ctrl = pyves.Controller()
        ctrl.readParameters("test/parameters.json")
        ctrl.prepareSimulation()
        ctrl.sample()
        lookup_table = ctrl.system.lookupTable

        r1 = ctrl.system.potentialEnergy(None)  
        r2 = ctrl.system.potentialEnergy(lookup_table)  
        r3 = ctrl.system.potentialEnergyConcurrent(None)    
        r4 = ctrl.system.potentialEnergyConcurrent(lookup_table)     



        print("\n\nINTERACTION RESULTS")
        print("Standard: Linear    ", r1)
        print("Lookup:   Linear    ", r2)
        print(f"Standard: {ctrl.system.threads} threads ", r3)
        print(f"Lookup:   {ctrl.system.threads} threads ", r4)


        # self.assertAlmostEqual()



    def test_interaction_benchmark(self):
        
        import pyves
        ctrl = pyves.Controller()
        ctrl.readParameters("test/parameters.json")

        setup = '''
import pyves
ctrl = pyves.Controller()
ctrl.readParameters("test/parameters.json")
ctrl.prepareSimulation()
ctrl.sample()
lookup_table = ctrl.system.lookupTable
'''
        standard = ''' 
a = ctrl.system.potentialEnergy(None)
'''

        lookup = ''' 
a = ctrl.system.potentialEnergy(lookup_table)
'''

        concurrent_standard = ''' 
a = ctrl.system.potentialEnergyConcurrent(None)
'''
        concurrent_lookup = ''' 
a = ctrl.system.potentialEnergyConcurrent(lookup_table)
'''

        num = 5000
        # timeit statement 
        import timeit
        print("\n\nINTERACTION TIMES")
        print("Standard: Linear    ", timeit.timeit(setup=setup, stmt=standard, number=num))
        print("Lookup:   Linear    ", timeit.timeit(setup=setup, stmt=lookup, number=num))

        print(f"Standard: {ctrl.system.threads} threads ", timeit.timeit(setup=setup, stmt=concurrent_standard, number=num))
        print(f"Lookup:   {ctrl.system.threads} threads ", timeit.timeit(setup=setup, stmt=concurrent_lookup, number=num))
        print("\n")
        print("\n")



if __name__ == '__main__':
    unittest.main()
