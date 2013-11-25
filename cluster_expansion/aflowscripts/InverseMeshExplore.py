import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice_vecs

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log
from numpy import matrix, transpose,rint,inner,multiply,size,argmin
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
    return abs(od)

def surfvol(vecs):
    vecs = transpose(vecs)
    u = norm(cross(vecs[0,:],vecs[1,:]))
    v = norm(cross(vecs[1,:],vecs[2,:]))
    w = norm(cross(vecs[2,:],vecs[0,:]))
    surfvol = 2*(u + v + w)/6/(abs(det(vecs)))**(2/3.0)
    return abs(surfvol)
    
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
def cosvecs(vec1,vec2):
    return dot(vec1,vec2)/norm(vec1)/norm(vec2)
        
def fillS(S,symops,nops):
    S2 = array(S)
    eps = 1.0e-6
    ''' Applies all symmetry operations, and chooses a new primitive vector that is 
    most orthogonal to the other(s)'''
    rotvecs = array_zeros((3,nops),dtype=np_float)
    dotvecs = array_zeros((nops),dtype=np_float)
#    if S2[0,1]==0 and S[1,1]==0 and S2[2,1]==0: 
    # we need to choose S2[:,1], the 2nd vector
    for iop in range(nops):
        rotvecs[:,iop] = transpose(dot(symops[:,:,iop],array(S2[:,0]))) # newvec = R S20; all 1-d arrays are horizontal
#        print 'rotvecs',iop
#        print rotvecs[:,iop]
        dotvecs[iop] = cosvecs(rotvecs[:,iop],S2[:,0])
#        print 'dotvecs',dotvecs[iop]
#    print '2nd vector',rotvecs[:,argmin(abs(dotvecs))]
    op2 = argmin(abs(dotvecs)); 
#    print 'Second vector found for op %i' % op2
    S2[:,1] =rotvecs[:,op2] #take the first most orthogonal vector
    #Find a 3rd vector
    done = False
    for iop in range(nops):
        rotvecs[:,iop] = transpose(dot(symops[:,:,iop],array(S2[:,1]))) # newvec = R S20; all 1-d arrays are horizontal
#        print 'rotvecs',iop
#        print rotvecs[:,iop]
        dot0 = cosvecs(rotvecs[:,iop],S2[:,0])
        dot1 = cosvecs(rotvecs[:,iop],S2[:,1])
#        print dot0,dot1
        if abs(abs(dot0)-1.0) > eps and abs(abs(dot1)-1.0) > eps: #have found independent vector
            S2[:,2] = rotvecs[:,iop]
#            print 'Third vector found for op %i' % iop; print S2
            done = True
            break
    if not done:
        print 'failed to find 3rd vector'  #one direction is uncoupled from the other two
    S = matrix(S2)
    return S    

def checksymmetry(latt,symops,nops):
    '''check that the lattice obeys all symmetry operations:  R.latt.inv(R) will give an integer matrix'''
    for iop in range(nops):
        lmat = array(latt)
        if det(lmat) == 0:
            return False
        mmat = dot(dot(inv(lmat),symops[:,:,iop]),lmat)
#        print 'mmat', iop
#        print mmat
        for i in range(3):
            for j in range(3):
                if abs(rint(mmat[i,j])-mmat[i,j])>1.0e-6:
#                    print iop, 'mmat[i,j]',mmat[i,j]                    
                    return False #jumps out of subroutine
    return True #passes check
    

##############################################################
########################## Script ############################

#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=200    


maxint = 8
#print 'Target N kpoints', Nmesh

M = zeros((3,3),dtype = np_int)
S = zeros((3,3),dtype = np_float)
B = lattice()
A = lattice()
K = lattice()

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

[symopsB,nopsB] = getGroup(B.vecs)
#print 'symmetry operations of B\n'
#for j in range(nopsB):
#    print j
#    op = matrix(symopsB[:,:,j])
#    print op
#find real lattice
A.vecs = trimSmall(inv(B.vecs).T)
A.det = det(A.vecs)
print 'A vectors';print A.vecs
print 'Det of A', A.det
print 'Orth Defect of A', orthdef(A.vecs)
print 'Surf/vol of A', surfvol(A.vecs)

[symopsA,nopsA] = getGroup(A.vecs)
if nopsA != nopsB: 
    print 'Number of operations different for A and B'
    sys.exit('stop')
else:
    nops = nopsA
#print 'symmetry operations of A\n'
#for k in range(nops):
#    print k
#    op = matrix(symopsA[:,:,k])
#    print op
    
#    sys.exit('stop')

#vary 
for a in range(maxint): # Diagonal elements of m can't be zero
    for b in range(maxint):
        for c in range(maxint):
            if a+b+c == 0:
                continue
            #create first S vector
            M[0,0]=a; M[0,1]=b; M[0,2]=c;
            S[:,0] = A.vecs*M.T[:,0]
#            print '1st vector';print S[:,0]
            #apply all symmetry operators, and find 2nd and 3rd vectors
            S = fillS(S,symopsA,nops)
            checksym = checksymmetry(S,symopsA,nops)
            if checksym: 
                K.vecs = transpose(inv(S))
                checksymB = checksymmetry(K.vecs,symopsB,nops)
                if checksymB: 
#                    print 'Obeys symmetry of lattice B:', checksymB 
                    print; print [a,b,c],'a,b,c'                              
                    print abs(det(S)/A.det),'Volume of S vs A:'
                    print round(surfvol(S),2),round(orthdef(S),2),'SV','OD' 
                            

sys.exit('stop')        
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
print 'Ended without minimum after maximum %i steps' % istep

    
