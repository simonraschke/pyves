from pyves import utility
import unittest
import pyves
import os
import pathlib
import numpy as np
import pandas as pd



class MainTest(unittest.TestCase):
    def test_simulation_run(self):
        print()
        try:
            os.remove("test/data.h5")
        except Exception as e:
            print("no data.h5 to remove",e)

        control = pyves.Controller()
        print("pyves version:", pyves.__version__)
        print("pyves concurrency model:", pyves.concurrency_model())
        
        control.readParameters("test/parameters.json")

        control.prepareSimulation()
        # self.assertEqual(len(control.system.particles), 100)
        box_dims = np.array([control.system.box.x, control.system.box.y, control.system.box.z])
        cells_per_dim = np.array(box_dims/control.cell_min_size).astype(int)
        self.assertEqual(len(control.system.cells), np.cumprod(cells_per_dim)[-1])

        # self.assertEqual(control.system.numParticlesInCells(), 100)
        self.assertTrue(control.system.assertIntegrity())

        control.sample(steps=500, timestats=True)

        for cell in control.system.cells:
            for particle in cell.particles:
                self.assertTrue(cell.contains(particle))

        control2 = pyves.Controller()
        control2.readParameters("test/parameters.json")
        control2.prepareSimulation()

        box_dims = np.array([control.system.box.x, control.system.box.y, control.system.box.z])
        cells_per_dim = np.array(box_dims/control.cell_min_size).astype(int)
        self.assertEqual(len(control2.system.cells), np.cumprod(cells_per_dim)[-1])

        control2.sample(timestats=True, analysis=True)
        control2.system.prepareSimulationStep()

        # self.assertEqual(control2.system.numParticlesInCells(), 100)
        for i,p in enumerate(control2.system.particles):
            self.assertTrue(p.assertIntegrity(), f"particle {i}")
        for i,c in enumerate(control2.system.cells):
            self.assertTrue(c.assertIntegrity(), f"cell {i}")
        
        self.assertTrue(control2.system.assertIntegrity())


        pyves.hdf2gro(
            inpath = os.path.join(control2.output["dir"], control2.output["filename"]),
            outpath = os.path.join(control2.output["dir"], "trajectory.gro"),
            atom_repr = dict(
                RAND = "O",
                S1 = "S",#
                S2 = "B",#
                P = "C"#
            ),
            prmspath = "test/parameters.json"
        )

        pyves.hdf2gro(
            inpath = os.path.join(control2.output["dir"], control2.output["filename"]),
            outpath = os.path.join(control2.output["dir"], "trajectory.gro"),
            atom_repr = dict(
                RAND = "O",
                S1 = "S",#
                S2 = "B",#
                P = "C"#
            ),
            box = control2.system.box,
            with_direction = True
        )
        self.assertTrue(os.path.exists(os.path.join(control2.output["dir"], "trajectory.gro")))
        
        pyves.writeVMDrc(outdir=control2.output["dir"], traj_file_name="trajectory.gro", num_bonds=len(control2.system.particles))
        self.assertTrue(os.path.exists(os.path.join(control2.output["dir"], ".vmdrc")))
        self.assertTrue(os.path.exists(os.path.join(control2.output["dir"], "vmd.rc")))
        
        pyves.analyzeTrajectory(
            # inpath = os.path.join(control2.input["dir"], control2.input["filename"]),
            # outpath = os.path.join(control2.output["dir"], "data2.h5"),
            prmspath = "test/parameters.json",
            timestats=True,
            threads=3
        )
        self.assertTrue(os.path.exists(os.path.join(control2.output["dir"], "data.h5")))
        # import pandas as pd
        # pd.options.display.max_rows = None
        # print(pd.read_hdf(os.path.join(control2.output["dir"], "data.h5"), key="/time1500")[["z","clustersize","surfacepot"]])
        # print(control2.system.interaction_surface_width)



    def test_complete_flow(self):
        try:
            os.remove("test/data.h5")
        except Exception as e:
            print("no data.h5 to remove",e)
        ctrl = pyves.Controller.Static("test/parameters.json")

        try:
            os.remove("test/states.db")
        except Exception as e:
            print("no states.db to remove",e)
        pyves.gatherStates("test", dbpath="test/states.db", filenames=["data.h5","gradient.h5"], threads=1, sys=True, clstr=True, key_prefix="/time", clstr_min_size=5)
        
        dfc = pyves.readStates("test/states.db", sql="SELECT * FROM cluster_states")
        self.assertIsInstance(dfc, pd.DataFrame)
        self.assertGreater(dfc.index.size, 0)
        self.assertIn("x_mean", dfc.columns)
        self.assertIn("size", dfc.columns)
        
        dfs = pyves.readStates("test/states.db", sql="SELECT * FROM system_states")
        self.assertIsInstance(dfs, pd.DataFrame)
        self.assertGreater(dfs.index.size, 0)
        self.assertIn("x", dfs.columns)
        self.assertIn("clustersize", dfs.columns)



    def test_system_gradient_flow(self):
        try:
            os.remove("test/gradient.h5")
        except Exception as e:
            print("no gradient.h5 to remove",e)



        times = np.linspace(0,100,10, dtype=int)
        attribute = "temperature"
        values = [0.2 for i in range(5)] + [0.3 for i in range(5)]
        
        ctrl = pyves.Controller.DynamicSystemFlow("test/gradient.json", attr=attribute, times=times, values=values)
        self.assertEqual(ctrl.time_actual, 100)
        
        for i, time in list(enumerate(times))[1:]:
            _, metadata = utility.h5load("test/gradient.h5", f"/time{time}")
            self.assertAlmostEqual(values[i-1], metadata["temperature"], 5)

        del ctrl


        times = np.linspace(101,200,11, dtype=int)
        attribute = "temperature"
        values = [0.2 for i in range(6)] + [0.3 for i in range(5)]
        
        ctrl = pyves.Controller.DynamicSystemFlow("test/gradient.json", attr=attribute, times=times, values=values)
        self.assertEqual(ctrl.time_actual, 200)
        
        for i, time in list(enumerate(times))[1:]:
            _, metadata = utility.h5load("test/gradient.h5", f"/time{time}")
            self.assertAlmostEqual(values[i-1], metadata["temperature"], 5)



if __name__ == '__main__':
    unittest.main()