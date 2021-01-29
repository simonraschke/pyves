import os
from os.path import dirname
import subprocess
import sys
import io
import re
import json
import numpy as np
import pandas as pd
from subprocess import check_output
from shutil import which#, copy2
from pathlib import Path    
import distutils.log
import distutils.dir_util
from .default_input import default_prms


def kwargs2string(**kwargs):
    if not kwargs:
        return ""
    else:
        return ", ".join(f'{x[0]}={x[1]!r}' for x in kwargs.items())[0:]


def slurmSubmitScript(
    filename = None,
    dirpath = None,
    mode = "new",
    prmspath = None,
    prms = None,
    partition = None,
    threads = None,
    memory = None,
    hours = None,
    hours_min = None,
    email = None,
    nice = None,
    group = None,
    requeue = None,
    signal_time = 1000,
    python_path = sys.executable,
    controller = "pyves.Controller.Static",
    controller_kwargs = dict(analysis=True),
    dryrun = False,
    mkdir = True,
    log = False
):

    """
    will return submit script path
    """
    
    # if isinstance(python_path, type(None)):
    #     python_path = sys.executable

    if hours_min == None:
        hours_min = hours

    assert(hours_min <= hours)

    distutils.log.set_verbosity(distutils.log.DEBUG)
    prmsfilename = "parameters.json"

    
    # make new directory
    dirpath = os.path.realpath(dirpath)
    if log: print("trying dir:", dirpath)
    if os.path.exists(dirpath):
        dirpath += "_0"
        if log: print("exists. trying dir:", dirpath)
    _it = 1
    while os.path.exists(dirpath):
        dirpath = dirpath.rsplit("_",1)[0] + f"_{_it}"
        if log: print("exists. trying dir:", dirpath)
        _it += 1
    if log: print("valid path:", dirpath)
    

    if mkdir and not dryrun:
        if log: print("mkdir:", dirpath)
        os.makedirs(dirpath)

    # verify that there is a directory
    if not mkdir and not dryrun:
        if not os.path.exists(dirpath):
            raise OSError(f"{dirpath} does not exist")



    # define target submit script name
    filepath = os.path.join(dirpath, filename)
    if log: print("trying filepath:", filepath)
    _it = 0
    while os.path.exists(filepath):
        filepath = os.path.join(dirpath, filename)
        filepath = filepath.rsplit(".",1)[0] + f"{_it}." + filepath.rsplit(".",1)[-1]
        if log: print("exists. trying file:", filepath)
        _it += 1
    if log: print("valid file:", filepath)
    

    if mode == "new":
        
        # make new parameters.json from dict
        if isinstance(prms, type({})):
            if log: print("prms is given. will ignore \"prmspath\"")
            prmspath = os.path.realpath(os.path.join(dirpath, prmsfilename))
                
            if Path(prmspath).is_file():
                if log: print(f"save the old parameter file {prmspath}")
                distutils.file_util.copy_file(prmspath, os.path.join(dirpath, "parameters_old.json"), verbose=1)
            
            setup_prms = default_prms
            setup_prms.update(prms)
            with open(prmspath, 'w') as fp:
                json.dump(setup_prms, fp=fp, indent=4)


        # copy existing *.json to parameters.json
        elif Path(os.path.realpath(prmspath)).is_file():
            prmspath = os.path.realpath(prmspath)
            
            if Path(prmspath).is_file():
                if log: print(f"save the old parameter file {os.path.join(dirpath, prmsfilename)}")
                if Path(os.path.join(dirpath, prmsfilename)).is_file():
                    distutils.file_util.copy_file(os.path.join(dirpath, prmsfilename), os.path.join(dirpath, "parameters_old.json"), verbose=1)
            
            if log: print(f"copy the input file")
            distutils.file_util.copy_file(prmspath, dirpath, verbose=1)

        else:
            raise IOError("neither prms nor prmspath were given")

    
    kwargs_string = kwargs2string(**controller_kwargs)
    command = f"{python_path} -c \\\"import pyves; {controller}('{prmsfilename}', {kwargs_string})\\\""
    days, hours = divmod(hours, 24)
    mindays, minhours = divmod(hours_min, 24)


    string = io.StringIO()
    print(f"#!{python_path}", file=string)
    print(f"", file=string)

    if group != None:
        print(f"#SBATCH -A {group}", file=string)
    
    if email != None:
        print(f"#SBATCH --mail-type=FAIL", file=string)
        print(f"#SBATCH --mail-user={email}", file=string)
    
    if requeue != None:
        print(f"#SBATCH --requeue", file=string)
        print(f"#SBATCH --open-mode=append", file=string)
    
    print(f"#SBATCH --export=NONE", file=string)
    print(f"#SBATCH --ntasks=1", file=string)
    print(f"#SBATCH --nodes=1", file=string)
    print(f"#SBATCH --cpus-per-task={threads}", file=string)
    print(f"#SBATCH --mem={memory}", file=string)
    print(f"#SBATCH --nice={nice}", file=string)
    print(f"#SBATCH --partition={partition}", file=string)
    print(f"#SBATCH --time={days:0>1}-{hours:0>2}:00:00", file=string)
    print(f"#SBATCH --time-min={mindays:0>1}-{minhours:0>2}:00:00", file=string)
    print(f"#SBATCH --signal=12@{int(signal_time)}", file=string)
    print(f"", file=string)
    
    print(f"import os, sys, signal", file=string)
    print(f"status = os.system(\"srun {command}\")", file=string)
    print(f"sys.stdout.flush()", file=string)
    print(f"sys.stderr.flush()", file=string)
    print(f"code = (status >> 8) & 0xFF", file=string)
    print(f"sig = status & 0xFF", file=string)
    print(f"if sig == 0 and code == signal.SIGUSR2.value:", file=string)
    print( "    rstatus = os.system(\"scontrol requeue {0}\".format(os.environ['SLURM_JOB_ID']))", file=string)
    print(f"    rcode = (rstatus >> 8) & 0xFF", file=string)
    print(f"    sys.exit(rcode)", file=string)
    
    if dryrun:
        print(string.getvalue())
    else:
        with open(filepath, "w") as slurmfile:
            print(string.getvalue(), file=slurmfile)
        
    return filepath



def sbatchSubmitScript(
    scriptpath,
    name = "pyves job",
    dryrun = False
):
    """
    will return job id
    """
    
    sbatchpath = which("sbatch")
    print("sbatchpath", sbatchpath)
    if sbatchpath == None:
        raise RuntimeError("sbatch not found")
    
    # start job from inside job dir
    cwd = os.path.realpath( os.getcwd() )
    if not os.path.exists(sbatchpath):
        raise RuntimeError(f"path does not exist: {sbatchpath}")
    try:
        if not dryrun:
            os.chdir(os.path.dirname( os.path.realpath( scriptpath )) )

            out = check_output([f"{sbatchpath}", "-J", f"{name}", f"{scriptpath}"])
            # out should look like
            # Submitted batch job 7154194

            os.chdir( cwd )

            try:
                id = int(re.search(r'\d+', out.decode("utf-8")).group())
            except Exception as e:
                raise RuntimeError("no jobid in sbatch output")
        else:
            id = 1337
    except Exception as e:
        os.chdir(cwd)
        raise e

    assert(id>15)
    return id



def sbatchGroTrajectory(prms_path, sbatch_kwargs, atom_repr=["C","O","S","H","B"], dryrun=False):
    sbatchpath = which("sbatch")
    # print("sbatchpath", sbatchpath)
    if sbatchpath == None:
        raise RuntimeError("sbatch not found")

    prms_path = os.path.realpath(prms_path)
    with open(prms_path, "r") as fp:
        prms = json.load(fp)

    dirpath = dirname(prms_path)
    datafile = os.path.join(dirpath, prms["output"].get("filename"))
    trajfile = os.path.join(dirpath, "trajectory.gro")

    particles = prms["system"]["particles"]
    names = list(particles.keys())
    representation = "dict("
    for i in range(len(names)):
        representation += f"{names[i]}='{atom_repr[i]}',"
    representation = representation[:-1]+")"
    
    sbatch_kwargs_string = " ".join([f"{arg}={val}" for arg, val in sbatch_kwargs.items()])
    cmd = f"{sbatchpath} {sbatch_kwargs_string} --wrap=\"{sys.executable} -c \\\"import pyves; pyves.hdf2gro(inpath='{datafile}', outpath='{trajfile}', atom_repr={representation}, prmspath='{prms_path}', with_direction='True')\\\"\""
    
    cwd = os.path.realpath( os.getcwd() )
    if not os.path.exists(sbatchpath):
        raise RuntimeError(f"path does not exist: {sbatchpath}")

    try:
        if not dryrun:
            os.chdir(dirpath)

            out = check_output(cmd, shell=True)
            # out should look like
            # Submitted batch job 7154194

            os.chdir( cwd )

            try:
                id = int(re.search(r'\d+', out.decode("utf-8")).group())
            except Exception as e:
                raise RuntimeError("no jobid in sbatch output")
        else:
            id = 1337
    except Exception as e:
        os.chdir(cwd)
        raise e

    assert(id>15)
    return id



def sbatchVMDRC(prms_path, sbatch_kwargs, dryrun=False, bond_radius=6):
    sbatchpath = which("sbatch")
    print("sbatchpath", sbatchpath)
    if sbatchpath == None:
        raise RuntimeError("sbatch not found")

    prms_path = os.path.realpath(prms_path)
    with open(prms_path, "r") as fp:
        prms = json.load(fp)

    cwd = os.path.realpath( os.getcwd() )
    if not os.path.exists(sbatchpath):
        raise RuntimeError(f"path does not exist: {sbatchpath}")

    dirpath = dirname(prms_path)
    datafile = os.path.join(dirpath, prms["output"].get("filename"))
    trajfile = os.path.join(dirpath, "trajectory.gro")
    
    sbatch_kwargs_string = " ".join([f"{arg}={val}" for arg, val in sbatch_kwargs.items()])
    cmd = f"{sbatchpath} {sbatch_kwargs_string} --wrap=\"{sys.executable} -c \\\"import pyves; import pandas as pd; num=pd.read_hdf('{datafile}', key='/time0').index.size; pyves.writeVMDrc(outdir='{dirpath}', traj_file_name='{trajfile}', num_bonds=num, bond_radius={bond_radius})\\\"\""

    try:
        if not dryrun:
            os.chdir(dirpath)

            out = check_output(cmd, shell=True)
            # out should look like
            # Submitted batch job 7154194

            os.chdir( cwd )

            try:
                id = int(re.search(r'\d+', out.decode("utf-8")).group())
            except Exception as e:
                raise RuntimeError("no jobid in sbatch output")
        else:
            id = 1337
    except Exception as e:
        os.chdir(cwd)
        raise e

    assert(id>15)
    return id