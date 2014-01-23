import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,nonzero,float64
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from itertools import combinations
#
#class min(object): #reciprocal lattice
#    def __init__(self):         
#        self.vecs = []
#        self.det = []
#        self.Nmesh = []
#        self.symops = []
#        self.nops = []

def changewhich(M,B,run):
    bestgrad = 0
    bestdel = zeros((3,3),dtype=int)
    Mold = M
    oldcost = cost(Mold,B,run)
#    print 'oldcost',oldcost
    bestindex=[-1,-1,0]#initialize
    for i in range(3):
        for j in range(3):
            M[i,j] += 1;delInc = cost(M,B,run)-oldcost; M[i,j] += -1
            if delInc < 0 and delInc < bestgrad: bestindex = [i,j,1];bestgrad = delInc
            M[i,j] += -1;delDec = cost(M,B,run)-oldcost;M[i,j] += 1;
            if delDec < 0 and delDec < bestgrad: bestindex = [i,j,-1];bestgrad = delDec
#            print i,j, delInc, delDec
    return bestindex

def findmin(M,B,run):
    '''Finds minimum cost for the lattice by varying the integers of m, in the element that gives steepest descent
    The 'run' indicates the cost function to use'''
    maxsteps = 1000
    istep = 1
    while istep<maxsteps:
        bestindex = changewhich(M,B,run)
        if bestindex[2]==0:#found minimum
            newcost = cost(M,B,run)
            break
        else:
            M[bestindex[0],bestindex[1]] += bestindex[2]
            newcost = cost(M,B,run)
        istep += 1
    if istep < maxsteps:
        print 'Found minimum after %i steps' % istep
#        print 'Best M'; print M
        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));            
    else:
        print 'Ended without minimum after maximum %i steps' % istep
        sys.exit('Stop')
    return [M,K]

def cost(M,B,run):
    if isequal(det(M),0):
        sys.exit('Det(M) = 0 in cost function. Stop')
    else: 
        K = lattice()
        K.vecs = dot(B.vecs,inv(M));K.det = abs(det(K.vecs))
    if run == 0:
        K = lattice()
        K.vecs = dot(B.vecs,inv(M));K.det = abs(det(K.vecs))
        Nscale =1*.05; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        cost = surfvol(K.vecs)*(1+Ncost)
#        pf = packingFraction(K.vecs)
#        cost = (1/pf)*(1+Ncost)
    elif run == 1:
        Nscale =1*1.0; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        shapescale = 1* 0.5
        shapecost = shapescale * surfvol(K.vecs)
#        shapecost = shapescale * (1/packingFraction(K.vecs))
        cost = symmetryError(K.vecs,B)*(1+Ncost)*(1+shapecost)
    else:
       sys.exit('Error in cost function. Stop')      
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
    pfB = packingFraction(B.vecs)
    print 'Packing fraction of B:', round(pfB,4)
    
    [B.symops,B.nops] = getGroup(B.vecs)
    print 'Number of symmetry operations', B.nops
    print 'Lattice type:', latticeType(B.nops)
    #print 'symmetry operations of B\n'
    #for j in range(nopsB):
    #    print j
    #    op = array(symopsB[:,:,j])
    #    print op
    #find real lattice
#    A.vecs = transpose(inv(B.vecs))
#    A.det = det(A.vecs)
#    A.Nmesh = Nmesh
#    print;print 'A vectors';print A.vecs
#    print 'Det of A', A.det
#    print 'Orth Defect of A', orthdef(A.vecs)
#    print 'Surf/vol of A', surfvol(A.vecs)
#    print 'Packing fraction of A:', round(packingFractionBK.vecs),4)
#    
#    [A.symops,A.nops] = getGroup(A.vecs)
#    if A.nops != B.nops: 
#        sys.exit('Number of operations different for A and B; stop')
#    print 'symmetry operations R of A\n'
#    for k in range(A.nops):
#        print 
#        print k
#        op = array(A.symops[:,:,k])
#        print'symop R of A'; print trimSmall(op)
#        m = trimSmall(dot(dot(inv(A.vecs), A.symops[:,:,k]),A.vecs))          
#        print'symop m'; print m            

    print 'Best mesh ignoring symmetry'
    M = zeros((3,3),dtype=int)
    a = rint(Nmesh**(1/3.0)); c = a; f = int(Nmesh/a/c)
    M[0,0]= a
    M[1,1]= c
    M[2,2]= f
    print 'Starting M'
    print M    
    run = 0
    [M,K] = findmin(M,B,run)
    print M
    print 'Packing fraction:', round(packingFraction(K.vecs),4)

    print 'Near mesh with required symmetry'
    run = 1
    [M,K] = findmin(M,B,run)
    print M
    print 'Packing fraction:', round(packingFraction(K.vecs),4)    
 
#    print 'Near mesh with large packing fraction'
#    run = 0
#    [M,K] = findmin(M,B,run)   
#    print M
#    print 'Packing fraction:', round(packingFraction(K.vecs),4)
#        
#    print 'Near mesh with required symmetry'
#    run = 1
#    [M,K] = findmin(M,B,run)

    print 'Symmetry of final mesh:',checksymmetry(K.vecs,B)
    if checksymmetry(K.vecs,B):
        print M
        print 'K mesh'; print K.vecs
        K.det = abs(det(K.vecs))
        print 'N of mesh', B.det/K.det, 'vs target', B.Nmesh
        print round(surfvol(K.vecs),4),round(orthdef(K.vecs),4),'SV of mesh,','OD'
        pf1 = packingFraction(K.vecs); pfmax = pf1
        print 'Packing fraction', round(pf1,4)
        print 'Lattice type:', latticeType(B.nops)
        
        print; print 'Try FCC-like substitution'
        kmesh2 = zeros((3,3),dtype = float)
        scale = 2/4**(1/3.0)
        kmesh2[:,0] = K.vecs[:,1]/scale + K.vecs[:,2]/scale
        kmesh2[:,1] = K.vecs[:,2]/scale + K.vecs[:,0]/scale
        kmesh2[:,2] = K.vecs[:,0]/scale + K.vecs[:,1]/scale 
        print kmesh2
        if checksymmetry(kmesh2,B):
            print 'N of mesh', B.det/det(kmesh2)
            SV = surfvol(kmesh2)
            print round(surfvol(kmesh2),4),round(orthdef(kmesh2),4),'SV of mesh,','OD'
            pf = packingFraction(kmesh2)
            print 'Packing fraction', round(pf,4)
            if pf > pfmax:
                K.vecs = kmesh2
                print 'FCC-like substitution is better'
                pfmax = pf

#        print; print 'Try BCC-like substitution'
#        kmesh2 = zeros((3,3),dtype = float)
#        scale = 2/2**(1/3.0)
#        kmesh2[:,0] = K.vecs[:,1]/scale + K.vecs[:,2]/scale - K.vecs[:,0]/scale
#        kmesh2[:,1] = K.vecs[:,2]/scale + K.vecs[:,0]/scale - K.vecs[:,1]/scale
#        kmesh2[:,2] = K.vecs[:,0]/scale + K.vecs[:,1]/scale - K.vecs[:,2]/scale
#        print kmesh2
#        if checksymmetry(kmesh2,B):
#            print 'N of mesh', B.det/det(kmesh2)
#            SV = surfvol(kmesh2)
#            print round(surfvol(kmesh2),4),round(orthdef(kmesh2),4),'SV of mesh,','OD'
#            pf = packingFraction(kmesh2)
#            print 'Packing fraction', round(pf,4)
#            if pf > pfmax:
#                K.vecs = kmesh2
#                print 'BCC-like substitution is better'
#                pfmax = pf                
                
        print
        print 'Final K mesh'; print K.vecs
        print 'Packing fraction', round(pfmax,4), 'vs original B', round(pfB,4)       
                
                    
        
        
        
    else:
        print'K mesh fails symmetry'  
        sys.exit('Stop')
