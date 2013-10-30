#!/usr/bin/python
''' Assumes vasp is ready to do mink reduction and chooses mesh. This overwrites puts "Minkowski Monkhorst-Pack" in place of "Monkhorst-Pack" 
in KPOINTS.  Replaces nx, ny, nz with "number of kpts per reciprocal atom"
'''
import sys,os,subprocess,time
import numpy as np
from numpy import pi
import kmeshroutines as km


km.checkq('bch') #loops for days
