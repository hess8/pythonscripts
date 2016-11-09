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


def joinLists(toJoin):
    '''Joins lists in toJoin of the form [[sublist1],[sublist2],[sublist3]...].  List must have the
    same length, but can different length sublists'''    
    joinedList = [[]]*len(toJoin[0])
    for isublist in range(len(toJoin[0])):
        sublist=[]
        for ilist in range(len(toJoin)):
            sublist += toJoin[ilist][isublist]
        joinedList[isublist] = sublist
    return joinedList

def multiDelete(list_, args):
    indexes = sorted(args, reverse=True)
    for index in indexes:
        del list_[index]
    return list_

a = [[1],[2,3],[9]]
b = [[4,5],[],[10]]
c = [[7],[8],[11]]
z = joinLists([a,b,c])
print z 

z = multiDelete(z,[0])
print z