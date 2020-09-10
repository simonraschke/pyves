import os
import pandas as pd
import pyves
import re



def _verify_io_files(i, o, f):
    if not os.path.exists(i):
        raise OSError(f"Input path {i} does not exists.")

    if not os.path.exists(o):
        if f:
            print(f"Output path {o} exists, will overwrite.")
        else:
            raise OSError(f"Output path {o} exists, but not allowed to overwrite.")



def hdf2gro(inpath, outpath, atom_repr={}, box=[], time_range=[0,1e10], overwrite=True, group_prefix="/time", with_direction=False):
    input_path = os.path.abspath(inpath)
    output_path = os.path.abspath(outpath)
    
    _verify_io_files(input_path, output_path, overwrite)

    store = pd.HDFStore(input_path, mode="r")
    keys = sorted([s for s in store.keys() if s.startswith(group_prefix)], key=lambda x:int(re.findall(f"(?<={group_prefix})\d+", x)[0]))
    
    with open(output_path, "w") as GROFILE:
        snapshots = 0
        for key in keys:
            time = key.split(group_prefix)[-1]
            if min(time_range) <= int(time) <= max(time_range):
                df = store.get(key)
                # if df["name"].unique().size != len(atom_repr):
                #     names = df["name"].unique()
                #     raise Exception(f"df names {names} do not match the atom representations {atom_repr}")

                if with_direction:
                    raise NotImplementedError("Orienataion vector representation not impletented.")
                else:
                    counter = 1
                    print(f"t={time}", file=GROFILE)
                    print(f"{df.index.size:>5}", file=GROFILE)
                    for i,row in enumerate(df.itertuples()):
                        # print(row)
                        line = f"{i+1:>5}"
                        line += str(f"{row.name}         ")[:5]
                        line += str(f"{atom_repr[row.name]:>5}       ")[:5]
                        line += f"{counter:>5}"
                        line += f"{row.x:>8.3f}{row.y:>8.3f}{row.z:>8.3f}"
                        line += "  0.0000  0.0000  0.0000"
                        print(line, file=GROFILE)
                        counter += 1
                    if isinstance(box, pyves.BoxPBC) or isinstance(box, pyves.BoxNoPBC):
                        print(f"{box.x:.4f} {box.y:.4f} {box.z:.4f}" , file=GROFILE)
                    else:
                        print(f"{box[1]:.4f} {box[1]:.4f} {box[2]:.4f}" , file=GROFILE)