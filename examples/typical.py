import pyves
import os

print("pyves typical use example")
print("pyves version:", pyves.__version__)

print("create Controller")
ctrl = pyves.Controller()

print("read forcefield")
ctrl.readForceField(os.path.join(os.path.dirname(__file__),"typical_forcefield.json"))

print("read parameters")
ctrl.readParameters(os.path.join(os.path.dirname(__file__),"typical_parameters.json"))

print("prepare simulation")
ctrl.prepareSimulation()

print("start sampling")
ctrl.sample(timestats=True)

pyves.hdf2gro(
    inpath = os.path.join(ctrl.output["dir"], ctrl.output["filename"]),
    outpath = os.path.join(ctrl.output["dir"], "trajectory.gro"),
    atom_repr = dict(
        FRAME = "O",
        MOBILE = "S"#
    ),
    box = ctrl.system.box
)