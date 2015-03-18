'''see kmeshroutines...this is a copy of the symmetry portion for Erik'''


'''Convention here is COLUMNS of matrices as vectors'''
################# functions #######################
from numpy import array, cos, sin,arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import sign, matrix, transpose,rint,inner,multiply,size,argmin,argmax,round,ceil
from numpy import zeros, nonzero, float64, sort, argsort, mod, amin, amax
fprec=float64
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy import int as np_int
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int
import time, os, subprocess, sys

utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') 
#had to copy and rename Gus's routine to the one below because ctypes could never find the one with the right name
getLatticePointGroup = utilslib.symmetry_module_mp_get_pointgroup_

def getGroup(latt):
#    print "lattice in getGroup\n",latt
    N = 3*3*48
    opsOUT =(c_double * N)() 
    NopsOUT =c_int(0) 
    lattIN = load_ctypes_3x3_double(transpose(latt)) # for some reason going to Fortran gives the TRANSPOSE
    eps = 1.0e-4
    epsIN = c_double(eps)
    getLatticePointGroup(byref(lattIN),byref(opsOUT),byref(NopsOUT),byref(epsIN)) 
    nops = NopsOUT.value
    symopsB = trimSmall(unload_ctypes_3x3xN_double(opsOUT,nops))
    return [symopsB,nops]

def intsymops(A):
    '''finds integer symmetry operations in the basis of the vectors of lattice A, 
    given the cartesian operators A.symop'''
    mops = zeros((3,3,A.nops),dtype = float)
    for i in range(A.nops):
        mops[:,:,i] = dot(dot(inv(A.vecs),A.symops[:,:,i]),A.vecs)
    return trimSmall(mops)  

def readposcar(filename, path): 
    ''' Format is explicit lattice vectors, not a,b,c,alpha, beta, gamma. 
    Saves vectors as columns, not rows which POSCAR uses'''
    file1 = open(path+filename,'r')
    poscar = file1.readlines()
    poscar = nstrip(poscar)
#    print poscar
    file1.close()
    descriptor = poscar[0]
    scale = float(poscar[1].split()[0])
    reallatt = zeros((3,3),dtype=fprec)
    reallatt[:,0] = array(poscar[2].split()[0:3])
    reallatt[:,1] = array(poscar[3].split()[0:3])
    reallatt[:,2] = array(poscar[4].split()[0:3])
#    reallatt = reallatt.astype(fprec)
#    print reallatt
    if scale < 0:
        vol = det(reallatt)
        scale = (-scale/vol)**(1/3.0)  
    reallatt = scale*reallatt
#    print reallatt    
    scale = 1.0
    reciplatt = 2*pi*transpose(inv(reallatt))
    natoms = array(poscar[5].split(),dtype=int)
    totatoms=sum(natoms)
    positions = zeros((totatoms,3))
    postype = poscar[6].split()[0] #Direct or Cartesian
    whichatom = 0
    for natom in natoms:
        for i in range(natom):
            for k in [0,1,2]:
                positions[whichatom,k] = float(poscar[7+whichatom].split()[k])
            whichatom += 1
    totatoms=sum(natoms)
    return [descriptor, scale, reallatt, reciplatt, natoms, postype, positions]

def create_poscar(filename,descriptor, scale, latticevecs, natoms, type_pos, positions, path):
    '''Write lattice vectors to POSCAR as rows, contrary to our convention.  
    The positions vectors are in an Nx3 matrix, so as rows already'''
    poscar = open(path+filename,'w')
    poscar.write(descriptor+'\n')
    poscar.write(str(scale)+'\n')
    for i in [0,1,2]:
        poscar.write('%20.15f %20.15f %20.15f \n' % (latticevecs[0,i], latticevecs[1,i], latticevecs[2,i])) #this column becomes a row in POSCAR       
    for i in natoms:
        poscar.write(str(i)+'    ')
    poscar.write('\n')
    poscar.write(type_pos+'\n')
    where = 0
    for natom in natoms:
        for i in range(natom):
            poscar.write('%20.15f %20.15f %20.15f \n' % (positions[where,0],positions[where,1],positions[where,2]))
            where += 1
    poscar.close()
        
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
    a = zeros((3,3,nops),dtype=fprec)
    ielement = 0
    for i in range(3):
        for j in range(3):
            for k in range(nops):                          
                a[i,j,k] = OUT[ielement]
                ielement += 1                 
    return a

def mink_reduce(a,eps):
    """Reduce the basis to the most orthogonal set.
       A Minkowski-reduced (via a "greedy algorithm basis) """
#        utilslib =  cdll.LoadLibrary('/Users/hart/codes/celib/trunk/libutils.so')
#    utilslib =  cdll.LoadLibrary('/fslhome/bch/cluster_expansion/theuncle/celib/trunk/libutils.so')
    utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so')
    ared =((c_double * 3) *3)()
#    mink = utilslib.vector_matrix_utilities_mp_minkowski_reduce_basis_ 
    mink = utilslib.vector_matrix_utilities_mp_minkowski_reduce_basis_     
    mink(byref(load_ctypes_3x3_double(a)),byref(ared),byref(c_double(eps)))
    ared2 = unload_ctypes_3x3_double(ared)   
    return ared2
