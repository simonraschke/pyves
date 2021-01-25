import os
from os.path import dirname
import sys
import io
import re
import json
import pyves
from subprocess import check_output
from shutil import which#, copy2
from pathlib import Path    
import distutils.log
import distutils.dir_util



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
            
            setup_prms = pyves.default_prms
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




    # if mode == "new":
    #     if log: print("mode new")
    #     # assert (not isinstance(prms, type(None))) or (not isinstance(prmspath, type(None)))
        
    #     # if (not isinstance(prms, type(None))) and (not isinstance(prmspath, type(None))):
    #     #     raise AttributeError(f"must give either prms or prmspath in mode new {prmspath}")
    #     # if (isinstance(prms, type(None))) and (isinstance(prmspath, type(None))):
    #     #     raise AttributeError(f"got prms and prmspath in mode new. please choose only one")
        
    #     # save the old parameter file
    #     if Path(os.path.join(dirpath, prmsfilename)).is_file():
    #         if log: print(f"save the old parameter file {os.path.join(dirpath, prmsfilename)}")
    #         distutils.file_util.copy_file(os.path.join(dirpath, prmsfilename), os.path.join(dirpath, "parameters_old.json"), verbose=1)

    #     # copy parameters file if path is given and exists
    #     if (not isinstance(prmspath, type(None))) and os.path.exists(prmspath):
    #         if log: print(f"copy the input file")
    #         distutils.file_util.copy_file(prmspath, dirpath, verbose=1)
        
    #     # when a prms dict is given, make input file from that
    #     elif isinstance(prms, type({})):
    #         default_prms = pyves.default_prms
    #         default_prms.update(prms)
    #         with open(os.path.join(dirpath, 'parameters.json'), 'w') as fp:
    #             json.dump(default_prms, fp)



    # elif mode == "restart":
    # elif mode == "gro":e
    


    # if os.path.dirname(os.path.realpath(prmspath)) != dirpath and not dryrun:
    #     if not os.path.isdir(dirpath):
    #         raise OSError(f"cant copy {filename} to {dirpath} . dir does not exist")
    #     if log: print("copy", prmspath, "to", dirpath)
    #     if isinstance(prms, type({})):
            
    #     elif update_prms_file:
    #         distutils.file_util.copy_file(prmspath, dirpath, update=1, verbose=1)
    #     # elif not os.path.exists(os.path.join(dirpath, filename)):
    #     #     distutils.file_util.copy_file(prmspath, dirpath, verbose=1)

    
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
    print(f"#SBATCH --signal=12@300", file=string)
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



def sbatchCommand(cmd, dirpath, sbatch_kwargs, dryrun=False):
    """
    will return job id
    """
    
    sbatchpath = which("sbatch")
    print("sbatchpath", sbatchpath)
    # if sbatchpath == None:
    #     raise RuntimeError("sbatch not found")

    sbatch_string = " ".join([f"{arg}={val}" for arg, val in sbatch_kwargs.items()])
    print(sbatch_string)
    print(f"{sbatchpath}", sbatch_string, cmd)
    
    cwd = os.path.realpath( os.getcwd() )
    if not os.path.exists(sbatchpath):
        raise RuntimeError(f"path does not exist: {sbatchpath}")
    try:
        if not dryrun:
            os.chdir(os.path.dirname( os.path.realpath( dirpath )) )

            out = check_output([f"{sbatchpath}", sbatch_string, cmd])
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