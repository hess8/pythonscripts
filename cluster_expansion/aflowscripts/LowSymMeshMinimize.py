'''Routines for minimizing S/V without regard to symmetry.  Useful for lowest symmetry crystals.'''

import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice,surfvol, orthdef

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, matrix, transpose,rint,inner,multiply
from numpy import zeros as array_zeros
from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv
from numpy import int as np_int
from numpy import float as np_float
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int

utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') 
#had to copy and rename Gus's routine to the one below because ctypes could never find the one with the right name


#

def Mfill(params):
#    print params
    [a,b] = params
    M = matrix((  #for BCT only
      [   a,  b,   -b],
      [   b,  a,   -b],
      [   0,   0,   a-b]
      ), dtype=np_int)  
    return M  


def changewhich(params,B):
    bestgrad = 0
    bestdel = zeros((3,3),dtype=np_int)
    Mold = Mfill(params)
    oldcost = cost(Mold,B)
    bestindex=[-1,0]#initialize (-1 not allowed, of course)
    params2 = params
    step = 1 
    for i,param in enumerate(params):         
            params2[i] += step;delInc = cost(Mfill(params2),B)-oldcost;params2[i] -= step           
            if delInc < 0 and delInc < bestgrad: bestindex = [i,1];bestgrad = delInc
            params2[i] -= step;delDec = cost(Mfill(params2),B)-oldcost;params2[i] += step
            if delDec < 0 and delDec < bestgrad: bestindex = [i,-1];bestgrad = delDec
#            print 'delInc, delDec',delInc, delDec
#    print bestgrad     
    return bestindex

def cost(M,B):
    if det(M) == 0:
        return 100
    else:
        K = lattice()
        K.vecs = B.vecs*M.I;K.det = det(K.vecs)
        if K.det == 0:
            return 100
        else:
            scale = .8
            Ncost = scale * abs((B.det/K.det -B.Nmesh))/B.Nmesh
        #    delN = 0.8  
        #    if B.det/K.det > (1+delN)*B.Nmesh or B.det/K.det < (1-delN)*B.Nmesh:
        #        cost = 1000 #reject this...number of mesh points too far off
        #    else:
        #        cost = surfvol(K.vecs)#; print'minimizing via surfac/volume'
#            cost = surfvol(K.vecs) + Ncost
            cost = surfvol(K.vecs)*(1+Ncost)        
            return(cost)

def findmin(B,Nmesh):
    print 'Target N kpoints', Nmesh
    a = int(rint(Nmesh**(1/3.0)));
    c = a
    b = 0
    params = [a,b]
    M = Mfill(params)
     
    print 'Starting M'
    print 'Det M', det(M)
    print M
    
#    
#    B.det = det(B.vecs)
#    B.Nmesh = Nmesh
#    print 'B vectors';print B.vecs
#    print 'Det of B', B.det
#    print 'Orth Defect of B', orthdef(B.vecs)
#    print 'Surf/vol of B', surfvol(B.vecs)
#    print 'symmetry operations of B\n'
#    print 'starting cost of M', cost(M,B)
#    
#    [symops,nops] = getGroup(B.vecs)
#    printmathematica = True 
#    
#    for k in range(nops):
#    #    print k
#        op = matrix(symops[:,:,k])
#    #    print op
#    #    print 'Inv(B) R B\n'
#    
#        if printmathematica:
#            print 'RMat[[All, All, %s]]=' % str(k+1)
#            print'{'
#            for i in range(3):
#                print '{%f,%f,%f}'% (op[i,0],op[i,1],op[i,2])
#                if i < 2: 
#                    print ','
#                else: 
#                    print '};'
#            
#            print 'mMat[[All, All, %s]]=' % str(k+1)
#            print'{'
#            for i in range(3):
#                print '{%f,%f,%f}'% (trimSmall(B.vecs.I * op * B.vecs)[i,0],trimSmall(B.vecs.I * op * B.vecs)[i,1],trimSmall(B.vecs.I * op * B.vecs)[i,2])
#                if i < 2: 
#                    print ','
#                else: 
#                    print '};'
#            print
#      
#            print
#      
#    #    sys.exit('stop')
    maxsteps = 10000
    istep = 1
    while istep<maxsteps:
        bestindex = changewhich(params,B)
        if bestindex[1]==0:#found minimum at previous M
            newcost = cost(M,B)        
            break
        else:
            params[bestindex[0]] += bestindex[1]
            M = Mfill(params)
            newcost = cost(M,B)
        istep += 1
    #    sys.exit('stop')        
    if istep < maxsteps:
        print
        print 'Found minimum after %i steps' % istep
        print 'An optimum M'; print M
        K = lattice();K.vecs = B.vecs*M.I; K.det = det(K.vecs)
        print 'An optimum K mesh\n', K.vecs  
        print 'Number of mesh points', B.det/K.det
        print 'Orth defect',orthdef(K.vecs)
        print 'Surface/vol', surfvol(K.vecs)    
        print 'Mesh vector lengths:'; print norm(K.vecs[:,0]),norm(K.vecs[:,1]),norm(K.vecs[:,2])
        k0 = K.vecs[:,0]; k1 = K.vecs[:,1]; k2 = K.vecs[:,2]
        cosgamma = k0.T*k1/norm(k0)/norm(k1)
        cosalpha = k1.T*k2/norm(k1)/norm(k2)
        cosbeta =  k2.T*k0/norm(k2)/norm(k0)       
        print 'Mesh vector cosines:'; print cosalpha, cosbeta, cosgamma
        print 'Check B = KM   \n', K.vecs*M  
        print '\n\n\nTranspose for use in VMD or POSCAR'
        print 'B'; print B.vecs.T
        print 'K'; print K.vecs.T
        print 
    else:
        print 'Ended without minimum after maximum %i steps' % istep

    
