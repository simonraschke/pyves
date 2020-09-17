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