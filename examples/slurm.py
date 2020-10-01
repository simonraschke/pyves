import pyves
import os, sys

path_to_submit_script = pyves.slurmSubmitScript(
    filename = "submit.py",
    dirpath = os.path.join(os.getcwd(), "slurmtest"),
    prmspath = os.path.join(os.getcwd(), "slurm.json"),
    partition = "sbatch partition",
    threads = 4,
    memory = "6G",
    hours = 1,
    hours_min = 1,
    email = "someone@somewhere.lol",
    nice = 1,
    requeue = True,
    python_path = sys.executable,
    controller = "pyves.Controller.StaticFlow",
    controller_kwargs = dict(analysis=True, analysis_inline=False),
    log = True
)

print("submit script path:", path_to_submit_script)

jobid = pyves.sbatchSubmitScript(path_to_submit_script, "slurm-example")
print("jobid:", jobid)