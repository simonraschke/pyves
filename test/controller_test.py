
import unittest
import pyves
import os
import numpy as np



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
            box = control2.system.box
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
        
        pyves.writeVMDrc(outdir=control2.output["dir"], traj_file_name="trajectory.gro", num_bonds=len(control2.system.particles))
        
        pyves.analyzeTrajectory(
            # inpath = os.path.join(control2.input["dir"], control2.input["filename"]),
            # outpath = os.path.join(control2.output["dir"], "data2.h5"),
            prmspath = "test/parameters.json",
            timestats=True,
            threads=3
        )



    def test_complete_flow(self):
        try:
            os.remove("test/data.h5")
        except Exception as e:
            print("no data.h5 to remove",e)
        ctrl = pyves.Controller.StaticFlow("test/parameters.json")



if __name__ == '__main__':
    unittest.main()