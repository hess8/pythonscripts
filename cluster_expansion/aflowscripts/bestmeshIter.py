import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction
    
from unconstrainedSVmin import unconstrainedSVsearch

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,nonzero,float64
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from itertools import combinations

def changewhich(M,B):
    bestgrad = 0
    bestdel = zeros((3,3),dtype=int)
    Mold = M
    oldcost = cost2(Mold,B)
#    print 'oldcost',oldcost
    bestindex=[-1,-1,0]#initialize
    for i in range(3):
        for j in range(3):
            M[i,j] += 1;delInc = cost2(M,B)-oldcost; M[i,j] += -1
            if delInc < 0 and delInc < bestgrad: bestindex = [i,j,1];bestgrad = delInc
            M[i,j] += -1;delDec = cost2(M,B)-oldcost;M[i,j] += 1;
            if delDec < 0 and delDec < bestgrad: bestindex = [i,j,-1];bestgrad = delDec
#            print i,j, delInc, delDec
    return bestindex

def cost2(M,B):
    if isequal(det(M),0):
        return 1000
    else:
        K = lattice()
        K.vecs = dot(B.vecs,inv(M));K.det = abs(det(K.vecs))
        Nscale =1*1.0; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        SVscale = 1* 1.2; SVcost = SVscale * surfvol(K.vecs)
        cost = symmetryError(K.vecs,B)*(1+Ncost)*(1+SVcost)
#        Nscale =1*.05; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
#        cost = surfvol(K.vecs)*(1+Ncost)        
        return(cost)  

def bestmeshIter(Blatt,Nmesh):
    '''
    Starts with MT made of eigenvectors of the m(R,A) operators. Explores noninteger changes in MT 
    to minimize the errors in symmetry and the cost in S/V and Nmesh 
    The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer 3x3 matrix
    So K = B Inv(M) .  Work in the inverse space of this problem, where we can work with M instead of Inv(M). 
    T(InvK) =  T(InvB)T(M).  
    
    Define S = T(InvK), and the real lattice A = T(InvB). So S = A T(M) is a superlattice on the real lattice.
           
    Minimization scheme'''
    
    ##############################################################
    ########################## Script ############################
   
    M = zeros((3,3),dtype = int)
    S = zeros((3,3),dtype = fprec)
    B = lattice()
    A = lattice()
    K = lattice()
       
    B.vecs = Blatt/2/pi  #Don't use 2pi constants in reciprocal lattice here.
    #############End BCT lattice
    eps = 1.0e-6
    B.det = det(B.vecs)
    B.Nmesh = Nmesh
    print 'B vectors';print B.vecs #
    #print 'B transpose'; print transpose(B.vecs)
    print 'Det of B', B.det
    print 'Orth Defect of B', orthdef(B.vecs)
    print 'Surf/vol of B', surfvol(B.vecs)
    
    [B.symops,B.nops] = getGroup(B.vecs)
    print 'Number of symmetry operations', B.nops
    print 'Lattice type:', latticeType(B.nops)
    #print 'symmetry operations of B\n'
    #for j in range(nopsB):
    #    print j
    #    op = array(symopsB[:,:,j])
    #    print op
    #find real lattice
    A.vecs = transpose(inv(B.vecs))
    A.det = det(A.vecs)
    A.Nmesh = Nmesh
    print;print 'A vectors';print A.vecs
    print 'Det of A', A.det
    print 'Orth Defect of A', orthdef(A.vecs)
    print 'Surf/vol of A', surfvol(A.vecs)
    
    [A.symops,A.nops] = getGroup(A.vecs)
    if A.nops != B.nops: 
        sys.exit('Number of operations different for A and B; stop')
#    print 'symmetry operations R of A\n'
#    for k in range(A.nops):
#        print 
#        print k
#        op = array(A.symops[:,:,k])
#        print'symop R of A'; print trimSmall(op)
#        m = trimSmall(dot(dot(inv(A.vecs), A.symops[:,:,k]),A.vecs))          
#        print'symop m'; print m            

    print 'First find min S/V ignoring symmetry'
    [M,K.vecs] = unconstrainedSVsearch(B)
    MT = trimSmall(dot(inv(K.vecs),B.vecs))
    print 'MT ignoring symmetry'; print MT
    print 'K.vecs';print K.vecs; 
#        for k in range(A.nops):

    maxsteps = 1000
    istep = 1
    while istep<maxsteps:
        bestindex = changewhich(M,B)
        if bestindex[2]==0:#found minimum
            newcost = cost2(M,B)
            break
        else:
            M[bestindex[0],bestindex[1]] += bestindex[2]
            newcost = cost2(M,B)
        istep += 1
    if istep < maxsteps:
        print
        print 'Found minimum after %i steps' % istep
        print 'Best M'; print M
        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));            
    else:
        print 'Ended without minimum after maximum %i steps' % istep
        sys.exit('Stop')
        
    print 'Symmetry of final mesh:',checksymmetry(K.vecs,B)
    if checksymmetry(K.vecs,B):
        print K.vecs
        K.det = abs(det(K.vecs))
        print 'N of mesh', B.det/K.det, 'vs target', B.Nmesh
        SV = surfvol(K.vecs)
        print round(surfvol(K.vecs),4),round(orthdef(K.vecs),4),'SV of Q2,','OD'
        pf = packingFraction(K.vecs)
        print 'packing fraction', pf
        print 'Lattice type:', latticeType(B.nops)
    else:
        print'K mesh fails symmetry'  
        sys.exit('Stop')
