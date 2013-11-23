import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice_vecs

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log
from numpy import matrix, transpose,rint,inner,multiply, size
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

def unload_ctypes_3x3xN_double(OUT,nops):
    """Take a ctypes 1-dim array and load it into a 3x3xnops python list.  Couldn't get 
    the 3x3xN to work"""
    a = array_zeros((3,3,nops),dtype=np_float)
    ielement = 0
    for i in range(3):
        for j in range(3):
            for k in range(nops):                          
                a[i,j,k] = OUT[ielement]
                ielement += 1                 
    return a

def getGroup(latt):
#    print "lattice in getGroup\n",latt
    N = 3*3*48
    opsOUT =(c_double * N)() 
    NopsOUT =c_int(0) 
    lattIN = load_ctypes_3x3_double(latt.T) # for some reason going to Fortran gives the TRANSPOSE
    eps = 1.0e-4
    epsIN = c_double(eps)
    getLatticePointGroup(byref(lattIN),byref(opsOUT),byref(NopsOUT),byref(epsIN)) 
    nops = NopsOUT.value
    symopsB = unload_ctypes_3x3xN_double(opsOUT,nops)
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
        
def newvec(S,symops,nops):
    ''' Applies all symmetry operations, and chooses a new primitive vector that is 
    most orthogonal to the other(s)'''
    rotvecs = array_zeros((3,nops),dtype=np_float)
    print  rotvecs[:,1]
    print  symops[:,1,1]
    if S[0,1]==0 and S[1,1]==0 and S[2,1]==0: # we need to find S[:,1], the 2nd vector
        for iop in range(nops):
            rotvecs[:,iop] = transpose(dot(symops[:,:,iop],array(S[:,0]))) # newvec = R S0; all 1-d arrays are horizontal
            print 'rotvecs',iop
            print rotvecs[:,iop]
    return newvec
#    if S[:,2] = [0.0,0.0]

#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=200    

print 'Target N kpoints', Nmesh

M = zeros((3,3),dtype = np_int)
S = zeros((3,3),dtype = np_float)
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
for j in range(nopsB):
    print j
    op = matrix(symopsB[:,:,j])
    print op
#find real lattice
A.vecs = trimSmall(inv(B.vecs).T)
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
for a in range(1,maxint): # Diagonal elements of m can't be zero
    for b in range(maxint):
        for c in range(maxint):
            #create first S vector
            M[0,0]=a; M[0,1]=b; M[0,2]=c;
            S[:,0] = A.vecs*M.T[:,0]
#            print 'S'; print S
            #apply all symmetry operators, and find a second vector
            S[:,1] = newvec(S,symopsA,nops)

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

    
