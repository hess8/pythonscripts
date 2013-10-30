from ctypes import byref, cdll, c_double
from numpy import array, zeros, transpose
    
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

def load_ctypes_3x3_double(IN):
    """Make a 3x3 array into the right thing for ctypes"""
    a = ((c_double * 3) *3)()
    for i in range(3):
        for j in range(3):
            a[i][j] = c_double(IN[i][j])
    return a

def unload_ctypes_3x3_double(OUT):
    """Take a ctypes array and load it into a 3x3 python list"""
    a = zeros((3,3))
    for i in range(3):
        for j in range(3):
            a[i][j] = OUT[i][j]
    return a

a = array(([0, .5, .5], [.5,0,.5],[.5,.5,0]), dtype=float)
a = transpose(array((  
  [   6.2133464115,   0.6213362303,  -0.6213362303],
  [   0.0000000000,   2.7083345463,   0.1425446910],
  [   0.0000000000,   0.0000000000,   2.7045807486]
  ), dtype=float))
eps = 1e-2
ared = array(mink_reduce(a,eps))
print 'Original lattice'
print a
print 'Reduced lattice'
print ared
