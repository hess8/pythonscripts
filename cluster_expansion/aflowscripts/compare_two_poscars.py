from ctypes import byref, cdll, c_double, c_bool,c_int,c_long
from numpy import array, zeros, transpose
from poscar import POSCAR as interpretposcar
import os, subprocess
    
def compare_structs(lat1,types1,pos1,natoms1,lat2,types2,pos2,natoms2,mapped,status,irot,skip_out,identical):
    """COMPARE ARBITRARY STRUCTURES
! This subroutine takes two structures, str1 and str2, and compares them to see if they are equivalent.
!  - Rescales the structures to have the same volume/atom  - Checks that both structures reside on the same
! underlying lattice, and that the atomic sites of str2 also lie on the lattice of str1
!  - Uses a set of spacegroup operations to find equivalent, but mis-oriented structures 
real(dp),intent(in):: LV1in(3,3), LV2in(3,3)    ! Lattice vectors for each structure
real(dp),intent(in):: aPos1in(:,:), aPos2in(:,:)! Atomic positions for each structure
integer,intent(in) :: aTyp1in(:), aTyp2in(:)    ! Atom types for each structure
real(dp), intent(in):: eps ! Finite precision tolerance
logical,intent(out):: mapped
logical,optional,intent(out):: identical
integer,optional,intent(out) :: status, irot, SMskipout
subroutine test_compare_arbitrary_structures(LV1in,aTyp1in,aPos1in,natoms1,LV2in,aTyp2in,aPos2in,natoms2,&
                                     eps,mapped,status,irot,SMskipout,identical)
                                     """
    utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') #note celib doesn't have this compare function in libutils.so
#    compare = utilslib.compare_structures_mp_compare_arbitrary_structures_ 
    test = utilslib.compare_structures_mp_test_compare_arbitrary_structures_     
#    ared =((c_double * 3) *3)()
    eps = 1.0e-4
#    compare_arbitrary_structures()

#    compare(byref(load_ctypes_3x3_double(lat1)),byref(load_ctypes_int_array(types1)),byref(load_ctypes_3x3xN_double(pos1)), \
#            byref(load_ctypes_3x3_double(lat2)),byref(load_ctypes_int_array(types2)),byref(load_ctypes_3x3xN_double(pos2)), \
#            byref(c_double(eps)),byref(c_bool(mapped)),byref(c_int(status)), \
#            byref(c_int(irot)),byref(c_int(skip_out)),byref(c_bool(identical)) )
    ntypes1 = len(types1)
    ntypes2 = len(types1)
    print
#    print lat1; print types1; print pos1
#    print load_ctypes_3x3_double(lat1); print load_ctypes_int_array(types1); print load_ctypes_3x3xN_double(pos1)
#    print 'len passed type1',len(load_ctypes_int_array(types1))
    test(byref(load_ctypes_3x3_double(lat1)),byref(load_ctypes_int_array(types1)),byref(c_int(ntypes1)), \
         byref(load_ctypes_3x3xN_double(pos1)), byref(c_int(natoms1)), \
            byref(load_ctypes_3x3_double(lat2)),byref(load_ctypes_int_array(types2)),byref(c_int(ntypes2)) \
            byref(load_ctypes_3x3xN_double(pos2)), byref(c_int(natoms2)), \
            byref(c_double(eps)),byref(c_bool(mapped)),byref(c_int(status)), \
            byref(c_int(irot)),byref(c_int(skip_out)),byref(c_bool(identical)) )    
    
    print mapped,status, irot, skip_out,identical 
    return 

def load_ctypes_3x3_double(IN):
    """Make a 3x3 array into the right thing for ctypes"""
    a = (c_double*3*3)()
    for i in range(3):
        for j in range(3):
            a[i][j] = c_double(IN[i][j])
    return a

def load_ctypes_3x3xN_double(IN):
    """Make a 3x3xN array into the right thing for ctypes"""
    N = IN.shape[1]
    a = (c_double *N*3)() #note different order N,3 vs standard matrix 3,N
#    print IN
#    print len(a)
#    print a
    for i in range(3):
        for j in range(N):
            a[i][j] = c_double(IN[i,j])            
    return a

def print_matrix(IN):
    N0 = IN.shape[0]
    N1 = IN.shape[1]
    for i in range(N0):
        for j in range(N1):
            print IN[i,j]
    
    
def load_ctypes_int_array(IN):
    """Make an integer array into the right thing for ctypes"""
    N = len(IN)
    a = (c_long * N)()
    for i in range(N):
        a[i] = c_long(IN[i])
        print 'int assigned',a[i]
    return a

def unload_ctypes_3x3_double(OUT):
    """Take a ctypes array and load it into a 3x3 python list"""
    a = zeros((3,3))
    for i in range(3):
        for j in range(3):
            a[i][j] = OUT[i][j]
    return a

def readposcar(filename): 
    file = open(filename,'r')
    rlines  = [i.strip() for i in file.readlines()] #assumes use of  < POSCAR
    return interpretposcar(lines=rlines)

###########################   Script  ############################
file1 = 'POSCAR1'
file2 = 'POSCAR2'
dir = '/fslhome/bch/cluster_expansion/alir/testf10/AlIr/f10/'
os.chdir(dir)
#read in both files
pos = readposcar(file1)
lat1 = array(pos.avecs)
types1 = array(pos.types)
natoms1 = len(pos.atoms)
pos1 = zeros((3,natoms1));pos2 = zeros((3,natoms1))
for i in range(natoms1):
    pos1[:,i]=pos.atoms[i].vector
print "Lattice 1:"; print lat1
print "Types 1:", types1
print "Positions 1:"; print pos1

pos = readposcar(file2)
lat2 = array(pos.avecs)
types2 = array(pos.types)
natoms2 = len(pos.atoms)
pos2 = zeros((3,natoms2))
for i in range(natoms2):
    pos2[:,i]=pos.atoms[i].vector
print "Lattice 2:"; print lat2
print "Types 2:", types2
print "Positions 2:"; print pos2

#compare the structures
mapped = False
identical = False
status = 0
irot = 0
skip_out = 0
compare_structs(lat1,types1,pos1,natoms1,lat2,types2,pos2,natoms2,mapped,status,irot,skip_out,identical)




#a = array(([0, .5, .5], [.5,0,.5],[.5,.5,0]), dtype=float)
#a = transpose(array((  
#  [   6.2133464115,   0.6213362303,  -0.6213362303],
#  [   0.0000000000,   2.7083345463,   0.1425446910],
#  [   0.0000000000,   0.0000000000,   2.7045807486]
#  ), dtype=float))
#eps = 1e-2
#ared = array(mink_reduce(a,eps))
#print 'Original lattice'
#print a
#print 'Reduced lattice'
#print ared
