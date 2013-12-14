'''Routines for minimizing S/V without regard to symmetry.  Useful for lowest symmetry crystals.'''

import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice,surfvol, orthdef

from kmeshroutines import lattice_vecs, lattice, surfvol, orthdef

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import matrix, transpose,rint,inner,multiply,size,argmin,nonzero,shape
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy import int as np_int
from numpy import float as np_float
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int

#from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, matrix, transpose,rint,inner,multiply
#from numpy import zeros as array_zeros
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
#from numpy.linalg import norm, det, inv
#from numpy import int as np_int
#from numpy import float as np_float
#from random import random, randrange
#from ctypes import byref, cdll, c_double, c_int

utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') 
#had to copy and rename Gus's routine to the one below because ctypes could never find the one with the right name

def searchmin(S,A):
    '''MT is transpose(M) '''
    MT = dot(inv(A.vecs),S)
    #determine how many empty rows there are in S:
    if norm(MT[:,0])==0: knownvecs = 0
    else:
        if norm(MT[:,1])==0: knownvecs = 1
        else:
            if norm(MT[:,2])==0: knownvecs = 2
    #Staring point of empty rows is a diagonal entry only
    for i in range(knownvecs,3):
        print A.Nmesh
        MT[i,i] = rint(A.Nmesh**(1/3.0))    
    #find the direction of steepest slope in the cost
    maxsteps = 10000
    istep = 1
    while istep<maxsteps:

        bestindex = changewhich(MT,knownvecs,A)
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
    
    
    return MT
                

def Mfill(params):
#    print params
    [a,b] = params
    M = matrix((  #for BCT only
      [   a,  b,   -b],
      [   b,  a,   -b],
      [   0,   0,   a-b]
      ), dtype=np_int)  
    return M  


def changewhich(MT,knownvecs,A):
    bestgrad = 0    
    bestdel = zeros((3,ncols),dtype=np_int)
    oldcost = cost(Mold,A)    
    bestchange=[[-1,-1],0]#initialize (-1 not allowed, of course)
    for i in range(3):
        for j in range(3-knownvecs):
            #find cost gradient for changes in integers by 1 
            MT[i,j] += 1;delInc = cost(MT,A)-oldcost; MT[i,j] -= 1           
            if delInc < 0 and delInc < bestgrad: bestindex = [[i,j],1];bestgrad = delInc
            MT[i,j] += -1;delDec = cost(MT,A)-oldcost; MT[i,j] -= -1           
            if delDec < 0 and delDec < bestgrad: bestindex = [[i,j],-1];bestgrad = delDec            
            print 'delInc, delDec',delInc, delDec            
    return bestindex

def cost(MT,A):
    if det(MT) == 0: sys.exit('Error det(M)=0')
    else:
        K = lattice()
        K.vecs = A.vecs*M.I;K.det = det(K.vecs)
        if K.det == 0:
            return 100
        else:
            scale = .8
            Ncost = scale * abs((A.det/K.det -A.Nmesh))/A.Nmesh
        #    delN = 0.8  
        #    if A.det/K.det > (1+delN)*A.Nmesh or A.det/K.det < (1-delN)*A.Nmesh:
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
