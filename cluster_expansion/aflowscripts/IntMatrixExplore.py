import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice_vecs

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, matrix, transpose,rint,inner
from numpy import zeros as array_zeros
from numpy.matlib import zeros #creates np.matrix rather than array
from numpy.linalg import norm, det
from numpy import int as np_int
from numpy import float as np_float
from random import random, randrange
from ctypes import byref, cdll, c_double

def get_class_members(klass):
    ret = dir(klass)
    if hasattr(klass,'__bases__'):
        for base in klass.__bases__:
            ret = ret + get_class_members(base)
    return ret


def uniq( seq ): 
    """ the 'set()' way ( use dict when there's no set ) """
    return list(set(seq))


def get_object_attrs( obj ):
    # code borrowed from the rlcompleter module ( see the code for Completer::attr_matches() )
    ret = dir( obj )
    ## if "__builtins__" in ret:
    ##    ret.remove("__builtins__")

    if hasattr( obj, '__class__'):
        ret.append('__class__')
        ret.extend( get_class_members(obj.__class__) )

        ret = uniq( ret )

    return ret

utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') 
#had to copy and rename Gus's routine to the one below because ctypes could never find the one with the right name
getLatticePointGroup = utilslib.symmetry_module_mp_get_pointgroup_

'''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer 3x3 matrix
So K = B Inv(M) 

         
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
def unload_ctypes_3x3x48_double(OUT):
    """Take a ctypes array and load it into a 3x3 python list"""
    a = array_zeros((3,3,48),dtype=np_float)

    for i in range(3):
        for j in range(3):
            for k in range(48):
                print [i,j,k], OUT[i][j][k]                
#                a[i,j,k] = OUT[i][j][k]
    return a

def getGroup(latt):
    opsOUT =(((c_double * 3) *3)*48)() 
    lattIN = load_ctypes_3x3_double(latt)
    eps = 1.0e-4
    epsIN = c_double(eps)
    getLatticePointGroup(byref(lattIN),byref(opsOUT),byref(epsIN))
    print opsOUT

#    for i in range(48):
#        print i
#        symops[:,:,i] = unload_ctypes_3x3_double(opsOUT[:][:][i])
    symops = unload_ctypes_3x3x48_double(opsOUT)
    return symops
    

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
    return bestindex

def cost(M,B):
    K = lattice()
    K.vecs = B.vecs*M.I;K.det = det(K.vecs)
    print 'symmetry operations of K\n', getGroup(K.vecs)
    if B.det/K.det > 1.1*B.Nmesh or B.det/K.det < 0.9*B.Nmesh:
        cost = 1000 #reject this...number of mesh points too far off
    else:
        cost = surfvol(K.vecs)#; print'minimizing via surfac/volume'
    return(cost)

def isinteger(x):
    return np.equal(np.mod(x, 1), 0)

class lattice(object): #reciprocal lattice
    def __init__(self):         
        self.vecs = []
        self.det = []
        self.Nmesh = []
#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=20
M = zeros((3,3),dtype=np_int)
print 'Target N kpoints', Nmesh
a = rint(Nmesh**(1/3.0)); c = a; f = int(Nmesh/a/c) 
M[0,0]= a
M[1,1]= c
M[2,2]= f
print 'Starting M'
print M
B =  lattice()

##############BCT lattice
alat = 2*sqrt(2)
clat = alat*4/3
B.vecs = matrix((  
  [   -alat/2,  alat/2,   alat/2],
  [   alat/2,  -alat/2,   alat/2],
  [   clat/2,   clat/2,   -clat/2]
  ), dtype=float)
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
#sys.exit('stop')
maxsteps = 10000
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
    sys.exit('stop')        
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

    
