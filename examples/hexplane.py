import pyves
import os, sys

ctrl = pyves.Controller.Static("hexplane.json", analysis=True)

print("make gro file from data")
pyves.hdf2gro(
    inpath = os.path.join(ctrl.output["dir"], ctrl.output["filename"]),
    outpath = os.path.join(ctrl.output["dir"], "trajectory.gro"),
    atom_repr = dict(
        CURVP = "O",
        PLANE = "C",
        CURVN = "N",
    ),
    box = ctrl.system.box,
    with_direction = True
)

print("write .vmdrc and vmd.rc")
pyves.writeVMDrc(outdir=ctrl.output["dir"], num_bonds=len(ctrl.system.particles))