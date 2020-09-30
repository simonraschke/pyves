import os
import sys
import io
import re
from subprocess import check_output
from shutil import which, copy2



def kwargs2string(**kwargs):
    if not kwargs:
        return ""
    else:
        return " ".join(f'{x[0]}={x[1]!r}' for x in kwargs.items())


def slurmSubmitScript(
    filename = None,
    dirpath = None,
    prmspath = None,
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
    controller = "pyves.Controller.StaticFlow",
    controller_kwargs = dict(analysis=True),
    dryrun = False,
    mkdir = True
):
    """
    will return submit script path
    """
    
    # if isinstance(python_path, type(None)):
    #     python_path = sys.executable

    if hours_min == None:
        hours_min = hours

    assert(hours_min <= hours)


    
    # make new directory
    dirpath = os.path.realpath(dirpath)
    if mkdir and not dryrun:
        if os.path.exists(dirpath):
            dirpath += "_0"
        _it = 1
        while os.path.exists(dirpath):
            dirpath = dirpath.split("_")[-2] + f"_{_it}"
            _it += 1
        os.makedirs(dirpath)

    # verify that there is a directory
    if not mkdir and not dryrun:
        if not os.path.exists(dirpath):
            raise OSError(f"{dirpath} does not exist")
    


    # define target submit script name
    filepath = os.path.join(dirpath, filename)
    _it = 0
    while os.path.exists(filepath):
        filepath = os.path.join(dirpath, filename)
        filepath = filepath.split(".")[-2]+f"{_it}."+filepath.split(".")[-1]
        _it += 1
    


    prmspath = os.path.realpath(prmspath)
    kwargs_string = kwargs2string(**controller_kwargs)
    command = f"{python_path} -c \\\"import pyves; {controller}('{prmspath}', {kwargs_string})\\\""
    days, hours = divmod(hours, 24)
    mindays, minhours = divmod(hours_min, 24)



    if os.path.dirname(os.path.realpath(prmspath)) != dirpath and not dryrun:
        if not os.path.isdir(dirpath):
            raise OSError(f"cant copy {filename} to {dirpath}  it does not exist")
        copy2(prmspath, dirpath)



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
    if not dryrun:
        os.chdir(os.path.dirname( os.path.realpath( scriptpath )) )

        out = check_output([f"{sbatchpath}", "-J", f"{name}", f"{scriptpath}"])
        # out should look like
        # Submitted batch job 7154194

        os.chdir( cwd )

        try:
            id = int(re.search(r'\d+', out).group())
        except Exception as e:
            raise RuntimeError("no jobid in sbatch output")
    else:
        id = 1337
        
    assert(id>15)
    return id