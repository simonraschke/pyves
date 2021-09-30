# pyves - python3 bindings for an easy use of vesicle2
# Copyright (C) 2020 Simon Raschke

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



import os
import pandas as pd
import pyves
import re
import json
from pathlib import Path



def _verify_io_files(i, o, f):
    if not os.path.exists(i):
        raise OSError(f"Input path {i} does not exists.")

    if not os.path.exists(o):
        if f:
            print(f"Output path {o} exists, will overwrite.")
        else:
            raise OSError(f"Output path {o} exists, but not allowed to overwrite.")



def hdf2gro(
        inpath, 
        outpath, 
        atom_repr={},
        uplus = None,
        uminus = None,
        box=None, 
        prmspath=None, 
        time_range=[0,1e10], 
        overwrite=True, 
        group_prefix="/time", 
        with_direction=True, 
        specific_key=None
    ):
    input_path = os.path.abspath(inpath)
    output_path = os.path.abspath(outpath)
    if not isinstance(prmspath, type(None)):
        prmspath = Path(os.path.abspath(prmspath))
        
    if (isinstance(box, type(None))):
        if prmspath.is_file():
            with open(prmspath, 'r') as prms_file:
                prms = json.loads(prms_file.read())
                box = pyves.BoxPBC([prms["system"]["box"]["x"], prms["system"]["box"]["y"], prms["system"]["box"]["z"]])
        else:
            raise FileNotFoundError(f"neither box nor prmspath \"{prmspath}\" are given")
    
    _verify_io_files(input_path, output_path, overwrite)

    store = pd.HDFStore(input_path, mode="r")
    keys = sorted([s for s in store.keys() if s.startswith(group_prefix)], key=lambda x:int(re.findall(f"(?<={group_prefix})\d+", x)[0]))
    if specific_key == "HEAD" or specific_key == "head":
        keys = [keys[-1]]
    elif specific_key != None:
        keys = [keys.index(specific_key)]
    
    with open(output_path, "w") as GROFILE:
        snapshots = 0
        for key in keys:
            time = key.split(group_prefix)[-1]
            if min(time_range) <= int(time) <= max(time_range):
                df = store.get(key).sort_values(by=['name'])

                atom = 1
                residue = 1
                if with_direction:
                    plus_vecs = df[["x", "y", "z"]].values + df[["ux", "uy", "uz"]].multiply(df["kappa"].div(2), axis="index").values
                    minus_vecs = df[["x", "y", "z"]].values - df[["ux", "uy", "uz"]].multiply(df["kappa"].div(2), axis="index").values
                    print(f"t={time}", file=GROFILE)
                    print(f"{df.index.size*2:>5}", file=GROFILE)
                    for i,row in enumerate(df.itertuples()):
                        plus = plus_vecs[i]
                        minus = minus_vecs[i]
                        # print(row)
                        line = f"{residue:>5}"
                        line += str(f"{row.name}         ")[:5]
                        if isinstance(uplus, type(None)):
                            line += str(f"{atom_repr[row.name]:>5}       ")[:5]
                        else:
                            line += str(f"{uplus:>5}       ")[:5]
                        line += f"{atom:>5}"
                        line += f"{plus[0]:>8.3f}{plus[1]:>8.3f}{plus[2]:>8.3f}"
                        line += "  0.0000  0.0000  0.0000"
                        print(line, file=GROFILE)
                        atom += 1

                        line = f"{i+1:>5}"
                        line += str(f"{row.name}         ")[:5]
                        if isinstance(uminus, type(None)):
                            line += str(f"{atom_repr[row.name]:>5}       ")[:5]
                        else:
                            line += str(f"{uminus:>5}       ")[:5]
                        line += f"{atom:>5}"
                        line += f"{minus[0]:>8.3f}{minus[1]:>8.3f}{minus[2]:>8.3f}"
                        line += "  0.0000  0.0000  0.0000"
                        print(line, file=GROFILE)
                        atom += 1
                        residue += 1
                else:
                    print(f"t={time}", file=GROFILE)
                    print(f"{df.index.size:>5}", file=GROFILE)
                    for i,row in enumerate(df.itertuples()):
                        # print(row)
                        line = f"{residue:>5}"
                        line += str(f"{row.name}         ")[:5]
                        line += str(f"{atom_repr[row.name]:>5}       ")[:5]
                        line += f"{atom:>5}"
                        line += f"{row.x:>8.3f}{row.y:>8.3f}{row.z:>8.3f}"
                        line += "  0.0000  0.0000  0.0000"
                        print(line, file=GROFILE)
                        atom += 1
                        residue += 1

                if isinstance(box, pyves.BoxPBC) or isinstance(box, pyves.BoxNoPBC):
                    print(f"{box.x:.4f} {box.y:.4f} {box.z:.4f}" , file=GROFILE)
                else:
                    print(f"{box[1]:.4f} {box[1]:.4f} {box[2]:.4f}" , file=GROFILE)



def writeVMDrc(outdir=".", traj_file_name="trajectory.gro", num_bonds=0, bond_radius=6):
    from pathlib import Path    
    
    vmdrc_paths = [
        os.path.abspath(os.path.join(outdir, "vmd.rc")),
        os.path.abspath(os.path.join(outdir, ".vmdrc"))
    ]
    for path in vmdrc_paths:
        with open(path, mode="w") as VMDRC:
            print("mol load gro",Path(traj_file_name).name, file=VMDRC)
            print("light 0 on", file=VMDRC)
            print("light 1 on", file=VMDRC)
            print("light 2 on", file=VMDRC)
            print("light 3 on", file=VMDRC)
            print("display nearclip set 0", file=VMDRC)
            print("axes location off", file=VMDRC)
            print("stage location off", file=VMDRC)
            print("menu main on", file=VMDRC)
            print("menu graphics on", file=VMDRC)
            print("display resize 1080 1080", file=VMDRC)
            print("display reposition", file=VMDRC)
            print("display projection perspective", file=VMDRC)
            # print("display rendermode GLSL", file=VMDRC)
            print("display cuedensity 0.12", file=VMDRC)
            print("color Display Background white", file=VMDRC)
            print("mol coloring 7 2 ResName", file=VMDRC)
            print(f"mol modstyle 0 0 Licorice {bond_radius} 10 30", file=VMDRC)
            print("mol modmaterial 0 0 AOChalky", file=VMDRC)
            print("# color definitions", file=VMDRC)
            print("color change rgb  0 0.07 0.20 0.48 ;# blue", file=VMDRC)
            print("color change rgb  1 0.70 0.20 0.10 ;# red", file=VMDRC)
            print("color change rgb  2 0.40 0.40 0.40 ;# gray", file=VMDRC)
            print("color change rgb  3 0.70 0.40 0.00 ;# orange", file=VMDRC)
            print("color change rgb  4 0.80 0.70 0.10 ;# yellow", file=VMDRC)
            print("color change rgb  7 0.13 0.47 0.04 ;# green", file=VMDRC)
            print("color change rgb  8 1.00 1.00 1.00 ;# white", file=VMDRC)
            print("color change rgb 10 0.10 0.70 0.80 ;# cyan", file=VMDRC)
            print("color change rgb 11 0.60 0.10 0.60 ;# purple", file=VMDRC)
            print("color change rgb 16 0.15 0.15 0.15 ;# black", file=VMDRC)
            print("after idle {", file=VMDRC)
            print("  pbc box -color black", file=VMDRC)
            print("  # set colors", file=VMDRC)
            print("  # create dummy molecule with one atom", file=VMDRC)
            print("  set mol [mol new atoms 1]", file=VMDRC)
            print("  set sel [atomselect $mol all]", file=VMDRC)
            print("  # add items to color categories", file=VMDRC)
            print("  # now we can define colors", file=VMDRC)
            # print("  color Name F 1", file=VMDRC)
            # print("  color Type F 1", file=VMDRC)
            print("  color Name S 23", file=VMDRC)
            print("  color Type S 23", file=VMDRC)
            print("  color Name C 16", file=VMDRC)
            print("  color Type C 16", file=VMDRC)
            print("  color Name O 1", file=VMDRC)
            print("  color Type O 1", file=VMDRC)
            print("  color Name N 7", file=VMDRC)
            print("  color Type N 7", file=VMDRC)
            print("  mol delete $mol", file=VMDRC)
            print("  # clean up", file=VMDRC)
            print("  $sel delete", file=VMDRC)
            print("#  set mol [mol new atoms 1]", file=VMDRC)
            print("#  set sel [atomselect $mol all]", file=VMDRC)
            print("  for {set x 0} {$x < ",num_bonds*2,"} {incr x} {", file=VMDRC)
            print("    set y [expr $x+1]", file=VMDRC)
            print("    set sel [atomselect top \"index $x $y\"]", file=VMDRC)
            print("    incr x 1", file=VMDRC)
            print("    set bonds [$sel getbonds]", file=VMDRC)
            print("    set ids [$sel get index]", file=VMDRC)
            print("    lassign $bonds atom1bonds atom2bonds", file=VMDRC)
            print("    lassign $ids atom1id atom2id", file=VMDRC)
            print("    lappend atom1bonds $atom2id", file=VMDRC)
            print("    lappend atom2bonds $atom1id", file=VMDRC)
            print("    $sel setbonds [list $atom1bonds $atom2bonds]" , file=VMDRC)
            print("    }", file=VMDRC)
            print("}", file=VMDRC)