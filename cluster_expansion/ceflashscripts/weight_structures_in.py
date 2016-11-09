import os, subprocess, sys, time 

from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,  \
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

def weightfunc(FE,FEmax,FEmin):
    maxratio = 1.0
    return str(1.0 + (maxratio-1)*(FEmax-FE)/(FEmax-FEmin))

'''Reads formation energies in structures.in.  If weighted, it changes the weights.  
If not, it adds a line for each weight'''

#dir = '/fslhome/bch/cluster_expansion/graphene/csv1-8W/'
dir = '/fslhome/bch/cluster_expansion/graphene/test1/'
infile = dir+'structures.in.441.FEsorted'
#infile = dir+'structures.in'
outfile = dir+'structures.in'
FEs = []
indices = []
lines = readfile(infile)
if 'noweights' in lines[1]: 
    weights_in = False
    lines[1] = 'weights\n'
else: 
    weights_in = True
#print 'weights_in', weights_in
#Find maximum and min formation energies
for i, line in enumerate(lines):
    if 'FE' in line:
        energy = float(line.split()[line.split().index('FE')+2].split(',')[0])  #FE = 1.2578,  etc.  remove any ","
        FEs.append(energy)
    if '#Energy' in line:
        indices.append(i+2) #where weights should appear 
FEmax = amax(FEs)
FEmin = amin(FEs)
indices = array(indices)

#Write weights
for i in range(len(FEs)):
    if weights_in:
        lines[indices[i]-2] = '#Energy and weight'+'\n'
        lines[indices[i]] = str(weightfunc(FEs[i],FEmax,FEmin))+'\n'
    else: #have to add lines
        lines.insert(indices[i],str(weightfunc(FEs[i],FEmax,FEmin))+'\n')
        indices += 1 #increment all indices by one to handle new lines
writefile(lines,outfile)

print 'Number of structures weighted', len(FEs)

print 'Done'
        



