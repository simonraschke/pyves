#!/usr/bin/python3
import pyves
import os

print("pyves typical use example")
print("pyves version:", pyves.__version__)

print("create Controller")
ctrl = pyves.Controller()

print("read forcefield")
ctrl.readForceField(os.path.join(os.path.dirname(__file__),"typical_forcefield.json"))

print("read parameters")
ctrl.readParameters(os.path.join(os.path.dirname(__file__),"typical_parameters.json"))

print("prepare simulation")
ctrl.prepareSimulation()

print("start sampling")
ctrl.sample()