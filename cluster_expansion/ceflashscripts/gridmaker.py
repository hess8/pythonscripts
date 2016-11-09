import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms, intsymops, create_poscar, searchsphere

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil,real,unravel_index, outer, fmod, amin, amax

from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations
from pylab import frange
#from ceScriptTools import runningJobs
from pylab import *
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy,vasputil_dosplot
 
dir = '/fslhome/bch/cluster_expansion/alal/cubic_al/mp_c1,3_nosymm/c1_4x4x4/'
#dir = '/fslhome/bch/cluster_expansion/alal/cubic_al/mp_c1,3_nosymm/c3_4x4x2/'
fold = 1
os.chdir(dir)

kmin = -0.5
kmax = 0.5
N = 2 # number in quadrant, not including gamma...must be even
dk = (kmax - kmin)/(2*N) #gamma centered

file1 = open('KPOINTS','w')
kpointsfile = []
kpointsfile.append('kpoints explicit grid BCH\n')
#kpointsfile.append('%i \n' % (N**3.0))
kpointsfile.append(str((2*N)**2 * 2*N/fold) +'\n' )
kpointsfile.append('Reciprocal \n')

count = 0
for i in range(-N+1,N+1):  #don't include edge states twice (-N+1), and add 1 to last one because of python's range features (N+1)
    for j in range(-N+1,N+1):
#        for l in range(N):
        for l in range(-N/fold+1,N/fold+1):
            print i,j,l
#            kpointsfile.append('%12.8f %12.8f %12.8f %12.8f\n' % (i*dk, j*dk, l*dk, 1.0))
            kpointsfile.append('%12.8f %12.8f %12.8f %12.8f\n' % (i*dk, j*dk, fold*l*dk, 1.0))
            count += 1
print count
kpointsfile.append('0.0 0.0 0.0\n' ) #shift
file1.writelines(kpointsfile) 
file1.close()
print 'done'



            
