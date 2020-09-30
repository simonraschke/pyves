import pyves
import os

print("pyves typical use example")
print("pyves version:", pyves.__version__)
print("pyves concurrency model:", pyves.concurrency_model())

print("create Controller")
ctrl = pyves.Controller()

print("read parameters")
ctrl.readParameters(os.path.join(os.path.dirname(__file__),"typical_parameters.json"))

print("prepare simulation")
ctrl.prepareSimulation()
ctrl.system.assertIntegrity()

print("start sampling")
ctrl.sample(timestats=True)

pyves.analyzeTrajectory(
    prmspath = os.path.join(os.path.dirname(__file__),"typical_parameters.json"),
    timestats = True
)
# ctrl = pyves.Controller.StaticFlow("typical_parameters.json")

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