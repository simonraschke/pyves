# pyves
## Python3 bindings for an easy use of vesicle2

### Installation

You may need to install
```bash
conda install gxx_linux-64
```

and afterwards 
```bash
python3 -m pip install --upgrade git+https://github.com/simonraschke/pyves.git
```
or
```bash
git clone https://github.com/simonraschke/pyves.git
cd pyves
git submodule update --init --recursive
pip install --upgrade .
```

### Usage
A simple example
```Python
import pyves

ctrl = pyves.Controller()
ctrl.readParameters("some_parameters.json") # see example in ./examples
ctrl.prepareSimulation()
ctrl.sample(timestats=True)
```

An even simpler example
```Python
ctrl = pyves.Controller.StaticFlow("some_parameters.json")
```

Using a control flow, that changes system parameters at certain time steps

```Python
ctrl = pyves.Controller.DynamicSystemFlow(
    "some_parameters.json", 
    attr = "temperature",     # the attribute to modify
    times = [0,10000,20000],  # change after these time values
    values = [0.25,0.24,0.23] # the values attr is set to
)
```


### Working with slurm

Direct control for single jobs.
All slurm jobs will be executed in the python environment you are currently in!
```Python
path_to_submit_script = pyves.slurmSubmitScript(
    filename = "submit.py",
    dirpath = "/path/to/my/job/dir/",
    prmspath = "/path/to/my/job/dir/parameters.json",
    partition = "sbatch partition",
    threads = 4,
    memory = "4G",
    hours = 24,
    hours_min = 23,                               # must be smaller than hours
    email = "someone@somewhere.lol",              # get an email if job failed
    nice = 0,
    requeue = True                                # will submit the job until it exits with not SIGUSR2 (12)
    # python_path = sys.executable,               # will use current environment executable by default
    # controller = "pyves.Controller.StaticFlow", # the default Controller
    # controller_kwargs = dict(analysis=True)     # controller argument dict
)

jobid = pyves.sbatchSubmitScript(path_to_submit_script, "some-job-name")
```

### Paramters
A parameter file `some_parameters.json` could look as follows (remember comments in json are not allowed and only included here for explanatory purposes).
```js
{
    "hardware" :                                // hardware parameters
    {   
        "threads" : 3                           // number of threads the program will spawn 
    },  
    "system" :                                  // Monte Carlo system parameters
    {   
        "temperature" : 0.20,                   // global system temperature
        "global_exchange_ratio" : 0.05,         // ratio of global exchange MC steps. 0.05 will generate an additional 5% global MC steps
        "global_exchange_epot_theshold": -1.0,  // epot threshold below which global exchange is possible (to avoid exchange of free beads)
        "box" :                                 // box size for periodic boundary distance calculations
        {   
            "x" : 20,   
            "y" : 15.5, 
            "z" : 18.3  
        },  
        "translation" :                         // translation of a particle per Monte Carlo step in Lennard Jones units
        {   
            "min" : 0.01,                       // minimum step width
            "max" : 0.2,                        // maximum step width
            "target" : 0.3,                     // target acceptance ratio
            "interval" : 1000                   // adjustment interval of the actual translations step width
        },  
        "rotation" :                            // rotation of a particle per Monte Carlo step in radiants
        {   
            "min" : 0.01,                       // minimum step width
            "max" : 0.2,                        // maximum step width
            "target" : 0.3,                     // target acceptance ratio
            "interval" : 1000                   // adjustment interval of the actual rotation step width
        },  
        "interaction" :                         // particle interaction parameters
        {   
            "cutoff" : 3.0                      // distance in Lennard Jones units above which the interaction energy is 0
        },  
        "particles" :                           // here, the actual particles get defined
        {   
            "A" :                               // an arbitrary name
            {   
                "number" : 1,                   // when dist==hexplane, this is the ratio of particles to similarly distributed
                "dist" : "hexplane",            // particles will form plane at box.z/2! 
                "sigma" : 1.0,                  // Lennard Joens parameter
                "kappa" : 1.0,                  // length of orientation vector (strength of anisotropy term)
                "epsilon" : 1.0,                // Lennard Joens parameter
                "gamma" : 0.1914626,            // optimum interaction angle in radiants
                "bound_translation" : null,     // restrict translation from its origin position to this
                "bound_rotation" : null         // restrict rotation from its origin orientation to this
            }   
            "B" :                               // an arbitrary name
            {   
                "number" : 1,                   // when dist==hexplane, this is the ratio of particles to similarly distributed
                "dist" : "hexplane",            // particles will form plane at box.z/2! 
                "sigma" : 1.0,                  // Lennard Joens parameter
                "kappa" : 1.0,                  // length of orientation vector (strength of anisotropy term)
                "epsilon" : 1.0,                // Lennard Joens parameter
                "gamma" : 0.1914626,            // optimum interaction angle in radiants
                "bound_translation" : null,     // restrict translation from its origin position to this
                "bound_rotation" : null         // restrict rotation from its origin orientation to this
            }   
            "S1" :                              // an arbitrary name
            {   
                "number" : 25,                  // the number of particles to be placed in the simulation box
                "dist" : "sphere200",           // a spherical distribution in the box center sized as big as an optimum 200 particle cluster
                "sigma" : 1.0,                  // Lennard Joens parameter
                "kappa" : 1.0,                  // length of orientation vector (strength of anisotropy term)
                "epsilon" : 1.0,                // Lennard Joens parameter
                "gamma" : 0.1914626,            // optimum interaction angle in radiants
                "bound_translation" : 0.1,      // restrict translation from its origin position to this
                "bound_rotation" : 0.3          // restrict rotation from its origin orientation to this
            },  
            "P" :                               // an arbitrary name
            {   
                "number" : 25,                  // the number of particles to be placed in the simulation box
                "dist" : "plane144",            // a planar distribution in the box center sized 12 x 12 = 144 in Lennard Jones units
                "sigma" : 1.0,                  // Lennard Joens parameter
                "kappa" : 1.0,                  // length of orientation vector (strength of anisotropy term)
                "epsilon" : 1.0,                // Lennard Joens parameter
                "gamma" : 0.1914626,            // optimum interaction angle in radiants
                "bound_translation" : 0.1,      // restrict translation from its origin position to this
                "bound_rotation" : 0.3          // restrict rotation from its origin orientation to this
            },  
            "RAND" :                            // an arbitrary name
            {   
                "number" : 175,                 // the number of particles to be placed in the simulation box
                "dist" : "random",              // a random particle distribtion
                "sigma" : 1.0,                  // Lennard Joens parameter
                "kappa" : 1.0,                  // length of orientation vector (strength of anisotropy term)
                "epsilon" : 1.0,                // Lennard Joens parameter
                "gamma" : 0.1914626,            // optimum interaction angle in radiants
                "bound_translation" : null,     // restrict translation from its origin position to this
                "bound_rotation" : null         // restrict rotation from its origin orientation to this
            }   
        }   
    },  
    "control" :                                 // algorithm control directly influencing performance
    {   
        "time_delta" : 1000,                    // time interval between output 
        "time_max" : 1000000,                   // maximum simulation time
        "cell_min_size" : 3.5,                  // minimum edge length per cell in Cell-Lists Algorithm
        "cell_update_interval" : 5,             // update interval of cell affiliation
        "neighbor_update_interval" : 5,         // update interval of neighbor lists
        "neighbor_cutoff" : 4.0                 // cutoff distance for neighbor list calculations
    },  
    "output" :                                  // 
    {   
        "dir" : "test",                         // directory
        "filename" : "data.h5",                 // filename with .h5
        "direct_analysis" : false               // analyze every simulation output directly (cannot be parallelized)
    },  
    "input" :   
    {   
        "dir" : "test",                         // directory
        "filename" : "data.h5",                 // filename with .h5
        "key" : "HEAD"                          // from where to restart simulation. HEAD is the maximum time stamp found
    },  
    "analysis": 
    {   
        "volume_max_grid" : 80,                 // volume calculation max grid per dimension (^3). less is faster and more is more accurate
        "volume_points_per_sigma" : 3,          // points per Lennard Jones sigma. less is faster and more is more accurate
        "DBSCAN_eps" : 1.2,                     // local density calulation cutoff. Is multiplied by the sigma of individual particles
        "neighbor_cutoff" : 1.2                 // neighbor calculation cutoff. Is multiplied by the sigma of individual particles
    }
}
```


### Citation
When using this software please cite

```
@article{raschke2019non,
    title={Non-equilibrium effects of micelle formation as studied by a minimum particle-based model},
    author={Raschke, Simon and Heuer, Andreas},
    journal={The Journal of chemical physics},
    volume={150},
    number={20},
    pages={204903},
    year={2019},
    publisher={AIP Publishing LLC}
}
```