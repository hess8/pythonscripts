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
from numpy import ceil,real,unravel_index, outer, fmod, amin, amax

#from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations
from pylab import frange


from ceroutines import readstructs,readunderlatt,scaleposcar, getscale, getline
''' 
struct_enum.in -> latt1, undetlying lattice
read struct list
for struct:
    create dir
    makestr.x -> poscar
    find J matrix by matrix inversion from latt1
    make a KPOINTS compatible with cubic and fcc mesh
'''

atomic = 'Al:Ir'
maindir = '/fslhome/bch/cluster_expansion/alir/enumtest/'
finaldir = '/fslhome/bch/cluster_expansion/alir/enumtest/structs/'
if not os.path.isdir(finaldir): os.system('mkdir %s' % finaldir)
structfile = '/fslhome/bch/cluster_expansion/alir/f1_5.dat'
enumfile = maindir + 'struct_enum.in'  
structchar = getline(0,enumfile).split('.')[-1][0]

#struct_enum.in -> latt1, underlying lattice A0
A0 = readunderlatt(enumfile)
print 'Underlying lattice';print A0
#enumerate lattice
os.system('enum.x %s' % enumfile)
#get scale from AFLOW
scale = getscale(atomic,structchar)
print scale
#read struct list
structs = readstructs(structfile)
#labeledstructs = [structchar + struct for struct in structs]
for struct in structs: #these are simply the numbers with no prefix
    dir = structchar + struct
    print dir
    if not os.path.isdir(finaldir+dir): os.system('mkdir %s' % finaldir+dir)
    os.chdir(finaldir+dir)
    os.system('makestr.x %s %s' % (maindir+'struct_enum.out', struct)) #creates poscar-like vasp.
    os.system('ls')
    os.system('mv vasp.* POSCAR')
    scaleposcar(scale)

print 'done'
