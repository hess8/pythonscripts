
import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms,points_in_ppiped #,three_perp,

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil,real,unravel_index

from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/') 
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  getms, writefermi, removezeros,getEf

'''reads chgcar from dir1 and multiplies by 2 (for superlattice of volume 2 along z axis), 
and writes it twice into CHGCAR for dir2 '''

dir1 = '/fslhome/bch/cluster_expansion/alal/cubic_al/mp_c1,3/c1_8x8x8/BANDS/'
dir2 = '/fslhome/bch/cluster_expansion/alal/cubic_al/mp_c1,3/c3_8x8x4/BANDS/'

lines1 = readfile(dir1 + 'CHGCAR')
lines2 = readfile(dir2 + 'CHGCAR')

print 'Done'