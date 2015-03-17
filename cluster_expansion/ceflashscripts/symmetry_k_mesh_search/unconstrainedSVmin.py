import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice, surfvol, orthdef, isequal,trimSmall
#from kmeshroutines import lattice,surfvol, orthdef

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, matrix, transpose,rint,inner,multiply
from numpy import equal, mod, logical_and
from numpy import zeros
from numpy.linalg import norm, det, inv
from numpy import int as np_int
from numpy import float as np_float
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int
def isinteger(x):
    return equal(mod(x, 1), 0)

def changewhich(M,B):
    bestgrad = 0
    bestdel = zeros((3,3),dtype=np_int)
    Mold = M
    oldcost = cost(Mold,B)
    bestindex=[-1,-1,0]#initialize
    for i in range(3):
        for j in range(3):
            M[i,j] += 1;delInc = cost(M,B)-oldcost; M[i,j] += -1
            if delInc < 0 and delInc < bestgrad: bestindex = [i,j,1];bestgrad = delInc
            M[i,j] += -1;delDec = cost(M,B)-oldcost;M[i,j] += 1;
            if delDec < 0 and delDec < bestgrad: bestindex = [i,j,-1];bestgrad = delDec
#            print i,j, delInc, delDec
    return bestindex

def cost(M,B):
    if isequal(det(M),0):
        return 100
#    if isequal(K.det,0):
#        return 100
    else:
        K = lattice()
        K.vecs = dot(B.vecs,inv(M));K.det = abs(det(K.vecs))
        Nscale =1*.05; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        cost = surfvol(K.vecs)*(1+Ncost)
        return(cost)     

    
def unconstrainedSVsearch(B):
    K = lattice();K.vecs = zeros((3,3),dtype=float)
    Nmesh = B.Nmesh
    M = zeros((3,3),dtype=np_int)
#    print 'Target N kpoints', Nmesh
    a = rint(Nmesh**(1/3.0)); c = a; f = int(Nmesh/a/c)
    M[0,0]= a
    M[1,1]= c
    M[2,2]= f
    print 'Starting M'
    print M
    maxsteps = 1000
    istep = 1
    while istep<maxsteps:
        bestindex = changewhich(M,B)
        if bestindex[2]==0:#found minimum
            newcost = cost(M,B)
            break
        else:
            M[bestindex[0],bestindex[1]] += bestindex[2]
            newcost = cost(M,B)
        istep += 1
    if istep < maxsteps:
        print
        print 'Found minimum after %i steps' % istep
#        print 'Best M'; print M
        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M))); 
       
#        K.det = det(K.vecs)
#        print 'Best K mesh\n', K.vecs
#        print 'Number of mesh points', B.det/K.det
#        print 'Minimum cost', newcost
#        print 'Orth defect',orthdef(K.vecs)
#        print 'Surface/vol', surfvol(K.vecs)
#        print 'Mesh vector lengths', norm(K.vecs[:,0]),norm(K.vecs[:,1]),norm(K.vecs[:,2])
    
    else:
        print 'Ended without minimum after maximum %i steps' % istep
    return [M,K.vecs]