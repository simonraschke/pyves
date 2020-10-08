import pyves
import os, sys

ctrl = pyves.Controller.StaticFlow("surface.json", analysis=False)

print("make gro file from data")
pyves.hdf2gro(
    inpath = os.path.join(ctrl.output["dir"], ctrl.output["filename"]),
    outpath = os.path.join(ctrl.output["dir"], "trajectory.gro"),
    atom_repr = dict(
        MOBILE = "S",
    ),
    box = ctrl.system.box,
    with_direction = True
)

print("write .vmdrc and vmd.rc")
pyves.writeVMDrc(outdir=ctrl.output["dir"], num_bonds=len(ctrl.system.particles))