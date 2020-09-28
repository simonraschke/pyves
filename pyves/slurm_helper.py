import os
import sys
import io



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
):
    assert(hours_min <= hours)

    # if isinstance(python_path, type(None)):
    #     python_path = sys.executable

    if hours_min == None:
        hours_min = hours

    filepath = os.path.join(os.path.abspath(dirpath), filename)
    _it = 0
    while os.path.exists(filepath):
        filepath = os.path.join(os.path.abspath(dirpath), filename)
        filepath = filepath.split(".")[-2]+f"{_it}."+filepath.split(".")[-1]
        _it += 1
    
    prmspath = os.path.abspath(prmspath)
    command = f"{python_path} -c \\\"import pyves; pyves.Controller.completeFlow('{prmspath}', analysis=False)\\\""
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
    
    with open(filepath, "w") as slurmfile:
        print(string.getvalue(), file=slurmfile)
        
    return filepath