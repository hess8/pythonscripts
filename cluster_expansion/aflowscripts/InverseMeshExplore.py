import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice_vecs

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
getLatticePointGroup = utilslib.symmetry_module_mp_get_pointgroup_

'''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer 3x3 matrix
So K = B Inv(M) .  Work in the inverse space of this problem, where we can work with M instead of Inv(M). 
T(InvK) =  T(InvB)T(M).  

Define S = T(InvK), and the real lattice A = T(InvB). So S = A T(M) is a superlattice on the real lattice.
       
Minimization scheme'''

def load_ctypes_3x3_double(IN):
    """Make a 3x3 array into the right thing for ctypes"""
    a = ((c_double * 3) *3)()
    for i in range(3):
        for j in range(3):
            a[i][j] = c_double(IN[i,j])
    return a

def unload_ctypes_3x3_double(OUT):
    """Take a ctypes array and load it into a 3x3 python list"""
    a = zeros((3,3))
    for i in range(3):
        for j in range(3):
            a[i][j] = OUT[i][j]
    return a
#def unload_ctypes_3x3x48_double(OUT):
#    """Take a ctypes array and load it into a 3x3 python list"""
#    a = array_zeros((3,3,48),dtype=np_float)
#    for i in range(3):
#        for j in range(3):
#            for k in range(48):
#                print [i,j,k], OUT[i][j][k]                
##                a[i,j,k] = OUT[i][j][k]
#    return a

def unload_ctypes_3x3x48_double(OUT,nops):
    """Take a ctypes 1-dim array and load it into a 3x3xnops python list"""
    a = array_zeros((3,3,nops),dtype=np_float)
    ielement = 0
    for i in range(3):
        for j in range(3):
            for k in range(nops):                          
                a[i,j,k] = OUT[ielement]
                ielement += 1                 
    return a

def getGroup(latt):
    print "lattice in getGroup\n",latt
    N = 3*3*48
    opsOUT =(c_double * N)() 
    NopsOUT =c_int(0) 
    lattIN = load_ctypes_3x3_double(latt.T) # for some reason going to Fortran gives the TRANSPOSE
    eps = 1.0e-4
    epsIN = c_double(eps)
    getLatticePointGroup(byref(lattIN),byref(opsOUT),byref(NopsOUT),byref(epsIN)) 
    nops = NopsOUT.value
    symopsB = unload_ctypes_3x3x48_double(opsOUT,nops)
    return [symopsB,nops]
    
def orthdef(latt): #columns as vectors
    od = norm(latt[:,0])*norm(latt[:,1])*norm(latt[:,2])/det(latt)
    return od

def surfvol(vecs):
    vecs = transpose(vecs)
    u = norm(cross(vecs[0,:],vecs[1,:]))
    v = norm(cross(vecs[1,:],vecs[2,:]))
    w = norm(cross(vecs[2,:],vecs[0,:]))
    surfvol = 2*(u + v + w)/6/(det(vecs))**(2/3.0)
    return surfvol
    
#def changewhich(M,B):
#    bestgrad = 0
#    bestdel = zeros((3,3),dtype=np_int)
#    Mold = M
#    oldcost = cost(Mold,B)
#    bestindex=[-1,-1,0]#initialize
#    for i in range(3):
#        for j in range(3):          
#            M[i,j] += 1;delInc = cost(M,B)-oldcost; M[i,j] += -1           
#            if delInc < 0 and delInc < bestgrad: bestindex = [i,j,1];bestgrad = delInc
#            M[i,j] += -1;delDec = cost(M,B)-oldcost;M[i,j] += 1;
#            if delDec < 0 and delDec < bestgrad: bestindex = [i,j,-1];bestgrad = delDec       
#    return bestindex

def Mfill(params):
#    print params
    [a,b] = params
    M = matrix((  #for BCT only
      [   a,  b,   -b],
      [   b,  a,   -b],
      [   0,   0,   a-b]
      ), dtype=np_int)  
    return M  

#def changewhich(M,B):
#    bestgrad = 0
#    bestdel = zeros((3,3),dtype=np_int)
#    Mold = M
#    oldcost = cost(Mold,B)
#    bestindex=[-1,-1,0]#initialize
#    for i in range(3):
#        for j in range(3):          
#            M[i,j] += 1;delInc = cost(M,B)-oldcost; M[i,j] += -1           
#            if delInc < 0 and delInc < bestgrad: bestindex = [i,j,1];bestgrad = delInc
#            M[i,j] += -1;delDec = cost(M,B)-oldcost;M[i,j] += 1;
#            if delDec < 0 and delDec < bestgrad: bestindex = [i,j,-1];bestgrad = delDec       
#    return bestindex

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

def isinteger(x):
    return np.equal(np.mod(x, 1), 0)

def trimSmall(mat):
    low_values_indices = abs(mat) < 1.0e-6
    mat[low_values_indices] = 0
    return mat

class lattice(object): #reciprocal lattice
    def __init__(self):         
        self.vecs = []
        self.det = []
        self.Nmesh = []

#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=200   #(a - b) (-b^2 + a c) for BCT 
#start b = 0 and c = a
#M = zeros((3,3),dtype=np_int)
print 'Target N kpoints', Nmesh
a = int(rint(Nmesh**(1/3.0)));
c = a
b = 0
params = [a,b]
M = Mfill(params)
 
print 'Starting M'
print 'Det M', det(M)
print M
B = lattice()
A = lattice()

##############BCT lattice
alat = 2*sqrt(2)
ca = 4/3.
clat = alat*ca
B.vecs = matrix((  
  [   -alat/2,  alat/2,   alat/2],
  [   alat/2,  -alat/2,   alat/2],
  [   clat/2,   clat/2,   -clat/2]
  ), dtype=float)
#print 'B vectors before inverse and transpose';print B.vecs
#B.vecs = trimSmall(inv(B.vecs)).T
#############End BCT lattice

############## Any lattice
#angle1 = 30
#crystal = [1,4,15,50,60,80] #trigonal [a,b,c,alpha,beta,gamma]
#B.vecs = transpose(lattice_vecs(crystal))
#############End BCT lattice
B.det = det(B.vecs)
B.Nmesh = Nmesh
print 'B vectors';print B.vecs
print 'Det of B', B.det
print 'Orth Defect of B', orthdef(B.vecs)
print 'Surf/vol of B', surfvol(B.vecs)
print 'symmetry operations of B\n'
[symopsB,nopsB] = getGroup(B.vecs)
for k in range(nopsB):
    print k
    op = matrix(symopsB[:,:,k])
    print op
#find real lattice
A.vecs = inv(B.vecs)
A.det = det(A.vecs)
print 'A vectors';print A.vecs
print 'Det of A', A.det
print 'Orth Defect of A', orthdef(A.vecs)
print 'Surf/vol of A', surfvol(A.vecs)
print 'symmetry operations of A\n'
[symopsA,nopsA] = getGroup(A.vecs)
if nopsA != nopsB: 
    print 'Number of operations different for A and B'
    sys.exit('stop')
else:
    nops = nopsA
for k in range(nops):
    print k
    op = matrix(symopsA[:,:,k])
    print op
    
#    sys.exit('stop')
maxint = 3
#vary 
for a in range(maxint):
    for b in range(maxint):
        for c in range(maxint):

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

    
