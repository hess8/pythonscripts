import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice_vecs, lattice, surfvol, orthdef
from LowSymMeshMinimize import findmin

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import matrix, transpose,rint,inner,multiply,size,argmin,nonzero
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
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
    a = zeros((3,3,nops),dtype=np_float)
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
    lattIN = load_ctypes_3x3_double(transpose(latt)) # for some reason going to Fortran gives the TRANSPOSE
    eps = 1.0e-4
    epsIN = c_double(eps)
    getLatticePointGroup(byref(lattIN),byref(opsOUT),byref(NopsOUT),byref(epsIN)) 
    nops = NopsOUT.value
    symopsB = unload_ctypes_3x3xN_double(opsOUT,nops)
    return [symopsB,nops]
    
def isinteger(x):
    return np.equal(np.mod(x, 1), 0)

def trimSmall(mat):
    low_values_indices = abs(mat) < 1.0e-6
    mat[low_values_indices] = 0.0
    return mat

def cosvecs(vec1,vec2):
    return dot(vec1,vec2)/norm(vec1)/norm(vec2)
        
def fillS(S,parentlatt):
    S2 = array(S)
    eps = 1.0e-6
    ''' Applies all symmetry operations, and chooses a new primitive vector that is 
    most orthogonal to the other(s)'''
    rotvecs = zeros((3,parentlatt.nops),dtype=np_float)
    dotvecs = zeros((parentlatt.nops),dtype=np_float)
    # we need to choose S2[:,1], the 2nd vector
    for iop in range(parentlatt.nops):
        rotvecs[:,iop] = transpose(dot(parentlatt.symops[:,:,iop],array(S2[:,0]))) # all 1-d arrays are horizontal
        #print 'rotvecs',iop
        #print rotvecs[:,iop]
        dotvecs[iop] = cosvecs(rotvecs[:,iop],S2[:,0])
        #print 'dotvecs',dotvecs[iop]
    if norm(abs(dotvecs)-1)<eps:# all rotated vectors are parallel or antiparellel to the first 
        return S2 #will have zero determinant, so skip this value of [a,b,c]
    #print '2nd vector',rotvecs[:,argmin(abs(dotvecs))]
    op2 = argmin(abs(dotvecs)); 
    #print 'Second vector found for op %i' % op2
    S2[:,1] =rotvecs[:,op2] #take the first most orthogonal vector
    #print 'S after 2nd vector'; print S2
    #Find a 3rd vector
    done = False
    for iop in range(parentlatt.nops):
        rotvecs[:,iop] = transpose(dot(parentlatt.symops[:,:,iop],array(S2[:,1]))) # newvec = R S20; all 1-d arrays are horizontal
        #print 'rotvecs',iop
        #print rotvecs[:,iop]
        dot0 = cosvecs(rotvecs[:,iop],S2[:,0])
        dot1 = cosvecs(rotvecs[:,iop],S2[:,1])
#        print dot0,dot1
        S2[:,2] = rotvecs[:,iop]
        if abs(det(S2)) >0.0 and abs(abs(dot0)-1.0) > eps and abs(abs(dot1)-1.0) > eps: #have found independent vector
            #print 'Third vector found for op %i' % iop; print S2
            done = True
            break
    if not done:
#        print 'failed to find 3rd vector via symmetry'  #one direction is uncoupled from the other two
        #choose a principal axis of the symmetry operators (rotation axis).      
        unitaxis = eigenvecfind(S2,parentlatt)
        #print 'unitaxis', unitaxis #we know the direction
        ghj = dot(inv(parentlatt.vecs),unitaxis) 
        #print 'ghj before normalization', ghj
        min1 = 1
        nonzero = abs(ghj) > 1.0e-6 #indices of nonzero
        ghj = ghj/min(abs(ghj[nonzero])); #should be array of integers in last row of M:  [g,h,j] for the smallest possible superlattice cell with vectors S2[:,0]), S2[:,1]) 
        #print 'inv(A)*unitaxis', ghj
#        sys.exit('stop')         
        #print dot(parentlatt.vecs,transpose(ghj))
        S2[:,2] = transpose(dot(parentlatt.vecs,transpose(ghj)))
        #print 'Third vector via eigenvector and minimizing S/V'; #print S2
        S2 = svmin_3rdvec(S2,parentlatt) #adjust the length of this last vector
        #print 'ending S in fillS'; print S2
        #
#        print S2

    S = trimSmall(S2)
    return S    

def svmin_3rdvec(S,parentlatt):
    '''Multiply third vector of S by integers and find the minimum S/Volume metric'''
    from copy import copy 
    sv = 100 #initialize too large
    trymax = 1000
    vec3min = copy(S[:,2])
    #print 'vec3min', vec3min
    for m in range(1,trymax):
        S[:,2] = float(m) * vec3min
        sv2 = surfvol(S)
        #print; print 'mult,sv',m,sv2
        #print vec3min
        #print'S'; print trimSmall(S)

        if sv2 > sv: #the last one was a minimum
            S[:,2] = (m-1) * vec3min
            #print 'best mult', m-1, S[:,2]
            return S
        else:
            sv = sv2
    sys.exit('Error in svmin_3rdvec; stop')   #sv kept growing  

def eigenvecfind(S,parentlatt):
    '''Find eigenvectors of symmetry operations (in cartesian basis).  A single eigenvector of eigenvalue 1 
    signifies a rotation operator. This returns a cartesian vector that corresponds to a rotation axis that is not in 
    the plane of the first two vectors of S'''
#    print 'S after array convert'; print S  
    eps = 1.0e-6
    
    for iop in range(parentlatt.nops): #loop for printing only
        evals,evecs = eig(parentlatt.symops[:,:,iop])     
        print; print iop; 
        print 'symop R'; print parentlatt.symops[:,:,iop]
#        print 'symop m'; print dot(dot(inv(parentlatt.vecs[:,:]),parentlatt.symops[:,:,iop]),parentlatt.vecs[:,:])    
#        print 'det(m)', det( dot(dot(inv(parentlatt.vecs[:,:]),parentlatt.symops[:,:,iop]),parentlatt.vecs[:,:]) )  
        print 'eigenvalues',evals
        print 'eigenvectors'; print evecs #loop for printing only
        
    for iop in range(parentlatt.nops):
        evals,evecs = eig(parentlatt.symops[:,:,iop])      
        nOne = 0
        for i, eval in enumerate(evals):
            if abs(eval - 1) < eps:
                nOne += 1
                indexOne = i
#                print 'indexOne',indexOne
        if nOne == 1:
#            print; print iop;print 'eigenvalues',evals
#            print 'eigenvectors'; print evecs 
            axis =  evecs[indexOne]
            unitaxis = axis/norm(axis) #unit vector
            #print "S";print S
            #print 'axis', axis          
            #test to see if it's independent of the first two axes
            #print 'cos angle between (cross(S[:,0],S[:,1]) and axis)', cosvecs(cross(S[:,0],S[:,1]),axis)
#            print 'cross(S[:,0],S[:,1]),axis', cross(S[:,0],S[:,1]),axis
            if abs(cosvecs(cross(S[:,0],S[:,1]),axis))>eps: 
                #print'found independent third vector'
                return unitaxis
#    #print 'Found no 3rd vector through eigenvectors of symmetry operators'
#    sys.exit('stop')
      
def checksymmetry(latt,parentlatt):
    '''check that the lattice obeys all symmetry operations of a parent lattice:  R.latt.inv(R) will give an integer matrix'''
    for iop in range(parentlatt.nops):
        lmat = array(latt)
        if det(lmat) == 0:
            return False
        mmat = dot(dot(inv(lmat),parentlatt.symops[:,:,iop]),lmat)
#        print 'mmat', iop
#        print mmat
        for i in range(3):
            for j in range(3):
                if abs(rint(mmat[i,j])-mmat[i,j])>1.0e-6:
#                    print iop, 'mmat[i,j]',mmat[i,j]                    
                    return False #jumps out of subroutine
    return True #passes test
    

##############################################################
########################## Script ############################

#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=200    

maxint = 2 
#print 'Target N kpoints', Nmesh

M = zeros((3,3),dtype = np_int)
S = zeros((3,3),dtype = np_float)
B = lattice()
A = lattice()
K = lattice()

##############BCT lattice
alat = 2*sqrt(5)
ca = 4/3.
clat = alat*ca
B.vecs = matrix((  
  [   -alat/2,  alat/2,   alat/2],
  [   alat/2,  -alat/2,   alat/2],
  [   clat/2,   clat/2,   -clat/2]
  ), dtype=float)


#B.vecs = matrix((  #C asis along x axis!####
#  [   -clat/2,  clat/2,   clat/2],
#  [   alat/2,  -alat/2,   alat/2],
#  [   alat/2,   alat/2,   -alat/2]
#  ), dtype=float)

#print 'B vectors before inverse and transpose';print B.vecs
#B.vecs = trimSmall(transpose(inv(B.vecs)))
#############End BCT lattice

############## Any lattice

#crystal = [1,1,sqrt(2),90,90,120] # [a,b,c,alpha,beta,gamma]
#crystal = [1,1,2,50,50,50] # [a,b,c,alpha,beta,gamma]
#B.vecs = transpose(lattice_vecs(crystal))
#############End BCT lattice
eps = 1.0e-6
B.det = det(B.vecs)
B.Nmesh = Nmesh
print 'B vectors';print B.vecs
print 'B transpose'; print transpose(B.vecs)
print 'Det of B', B.det
print 'Orth Defect of B', orthdef(B.vecs)
print 'Surf/vol of B', surfvol(B.vecs)

[B.symops,B.nops] = getGroup(B.vecs)
print 'Number of symmetry operations', B.nops
#print 'symmetry operations of B\n'
#for j in range(nopsB):
#    print j
#    op = matrix(symopsB[:,:,j])
#    print op
#find real lattice
A.vecs = trimSmall(transpose(inv(B.vecs)))
A.det = det(A.vecs)
print 'A vectors';print A.vecs
print 'Det of A', A.det
print 'Orth Defect of A', orthdef(A.vecs)
print 'Surf/vol of A', surfvol(A.vecs)

[A.symops,A.nops] = getGroup(A.vecs)
if A.nops != B.nops: 
    sys.exit('Number of operations different for A and B; stop')

print 'symmetry operations R of A\n'
for k in range(A.nops):
    print; print k
    op = matrix(A.symops[:,:,k])
    print trimSmall(op)
    m = trimSmall(dot(dot(inv(A.vecs[:,:]), A.symops[:,:,k]),A.vecs[:,:])  ) 
    [val,vecs]=eig(m)
    print 'symop m'; print m

    print 'det(m)', det(m)
    print 'eigen of m',val
    print vecs
    print vecs[:,0]/abs(vecs[:,0])[nonzero(vecs[:,0])].min()
    print vecs[:,1]/abs(vecs[:,1])[nonzero(vecs[:,1])].min()
    print vecs[:,2]/abs(vecs[:,2])[nonzero(vecs[:,2])].min()
    
    if k == 0:
        print vecs[0,2]*A.vecs[:,0]+vecs[1,2]*A.vecs[:,1]+vecs[2,2]*A.vecs[:,2]
                        
sys.exit('stop')
if A.nops < 4: #has only 2 symm operations. 
    findmin(B,Nmesh)
else:
    minSV = 1000 #initialize
    for a in range(maxint): # Diagonal elements of m can't be zero
        for b in range(maxint):
            for c in range(2*maxint):
                if a==0 and b==0 and c==0:
                    continue
                #create first S vector
                M[0,0]=a; M[0,1]=b; M[0,2]=c;
#                print; print [a,b,c],'a,b,c'
                S[:,0] = dot(A.vecs,transpose(M[:,0]))
                if norm(S[:,0])<eps:
                    continue # all zeros, so go to next a,b,c
                print '1st vector';print S[:,0]
                #apply all symmetry operators, and find 2nd and 3rd vectors
                S = fillS(S,A)
                if abs(det(S))<eps:
                    print 'S has zero det for ',[a,b,c]
                    continue
    #            print S
                checksym = checksymmetry(S,A)
                if checksym: 
                    K.vecs = transpose(inv(S))
                    checksymB = checksymmetry(K.vecs,B)
                    if checksymB: 
    #                    print 'Obeys symmetry of lattice B:', checksymB               
                        print; print [a,b,c],'a,b,c'
                        print 'S';print S
                        print abs(det(S)/A.det),'Volume of superlattice'
                        print round(surfvol(S),2),round(orthdef(S),2),'SV of superlattice','OD'  
                        print round(surfvol(K.vecs),2),round(orthdef(K.vecs),2),'SV of k-mesh','OD'  
                        print 'M matrix abc;def;ghj';print trimSmall(transpose(dot(inv(A.vecs),S)))
                        if surfvol(K.vecs)<minSV:
                            minSV = surfvol(K.vecs)
                            bestK = K.vecs
                            bestabc = [a,b,c]
    print'------------------------------------------------'
    print 'a,b,c',bestabc
    print 'Best k mesh'; print bestK
    print 'Best surface/volume:', round(minSV,2)

print 'Done'