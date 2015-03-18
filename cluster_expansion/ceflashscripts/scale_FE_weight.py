import os, subprocess, sys, time 

from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms,points_in_ppiped #,three_perp,

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,amax, amin, argmin,argmax,nonzero,float64, identity
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

def scaledFE(FE,FEmax,FEmin):
    scale = 5
    FEsuper = 1.1* FEmax
    return str(scale*(FE-FEsuper)) #scaled FEs should all be negative

'''Reads formation energies in structures.in or .holdout  
Ignores weights, and writes for energy a scaled formation energy that spreads differences in FE for low energy,
and compresses differences in FE for high FE'''

#dir = '/fslhome/bch/cluster_expansion/graphene/csv1-8W/'
dir = '/fslhome/bch/cluster_expansion/graphene/test1/'
#infile = dir+'structures.in.441.FEsorted'
infile = dir+'structures.holdout'
#infile = dir+'structures.in'
#outfile = dir+'structures.in'
outfile = dir+'structures.holdout'
FEs = []
indices = []
lines = readfile(infile)

#Find maximum and min formation energies
for i, line in enumerate(lines):
    if 'FE' in line:
        energy = float(line.split()[line.split().index('FE')+2].split(',')[0])  #FE = 1.2578,  etc.  remove any ","
        FEs.append(energy)
    if '#Energy' in line:
        indices.append(i) #where label 'Energy' appeared 
FEmax = amax(FEs)
FEmin = amin(FEs)
indices = array(indices)

#Write new Fes
for i in range(len(FEs)):
        lines[indices[i]] = '#Energy (scaled formation energy)'+'\n'
        lines[indices[i]+1] = str(scaledFE(FEs[i],FEmax,FEmin))+'\n'
writefile(lines,outfile)

print 'Number of structures with formation energies scaled', len(FEs)

print 'Done'
        



