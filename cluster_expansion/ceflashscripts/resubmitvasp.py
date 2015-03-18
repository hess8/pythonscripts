#!/usr/bin/python
import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms, intsymops, create_poscar, searchsphere

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil,real,unravel_index, outer, fmod, amin, amax, sign

#from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations
from pylab import frange


# ==================================================================================
# ==================================================================================
# ==================================================================================
#maindir = '/fslhome/bch/cluster_expansion/alir/enumtest/'
maindir = '/fslhome/bch/cluster_expansion/sisi/equivk/'

finalDir = maindir + 'structs.cubmesh/'
#finalDir = maindir + 'structs.fccmesh/'

testfilestr = 'slurm'
teststr = 'Disk quota exceeded'



os.chdir(finalDir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    print dir    
    os.chdir(finalDir+dir)
    files = os.listdir(os.getcwd())
    for file in files:
        if testfilestr in file:
            testfile = file    #go through the loop; will take the last one that has this string        
    file1 = open(testfile,'r')
    lines = file1.readlines()
    file1.close()
    for line in lines:
        if teststr in line:
            subprocess.call(['sbatch', 'vaspjob']) #!!!!!!! Submit jobs
#            print 'resubmitted', dir
            break
                
        

print 'done'
