import pyves
import os

pyves.Controller.printRuntimeInfo()

print("create Controller")
ctrl = pyves.Controller()

print("read parameters")
ctrl.readParameters(os.path.join(os.path.dirname(__file__),"typical.json"))

print("prepare simulation")
ctrl.prepareSimulation()
ctrl.system.assertIntegrity()

print("start sampling")
ctrl.sample(timestats=True)

pyves.analyzeTrajectory(
    prmspath = os.path.join(os.path.dirname(__file__),"typical.json"),
    timestats = True
)
# ctrl = pyves.Controller.StaticFlow("typical.json")

print("make gro file from data")
pyves.hdf2gro(
    inpath = os.path.join(ctrl.output["dir"], ctrl.output["filename"]),
    outpath = os.path.join(ctrl.output["dir"], "trajectory.gro"),
    atom_repr = dict(
        FRAME = "O",
        MOBILE = "S",
        MOBILE2 = "P",
    ),
    box = ctrl.system.box,
    with_direction = True
)

print("write .vmdrc and vmd.rc")
pyves.writeVMDrc(outdir=ctrl.output["dir"], num_bonds=len(ctrl.system.particles))