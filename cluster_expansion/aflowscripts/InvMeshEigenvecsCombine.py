import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, lattice_vecs, lattice, surfvol, orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs
from LowSymMeshMinimize import searchmin #,searchminRandomStart

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import matrix, transpose,rint,inner,multiply,size,argmin,nonzero
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy import int as np_int
from numpy import float as np_float
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int
from itertools import combinations

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
    

        
def findNextVec(S,parentlatt,which):
    ''' Applies all symmetry operations, and chooses a new primitive vector that is 
    most orthogonal to the other(s)'''    
    eps = 1.0e-6
    found = False
    newvec = array([0.,0.,0.])  
    rotvecs = zeros((3,parentlatt.nops),dtype=np_float)
    dotvecs = zeros((parentlatt.nops),dtype=np_float)
    dotvecs0 = zeros((parentlatt.nops),dtype=np_float)
    for iop in range(parentlatt.nops):
        rotvecs[:,iop] = transpose(dot(parentlatt.symops[:,:,iop],S[:,which-1])) # all 1-d arrays are horizontal
        dotvecs[iop] = cosvecs(rotvecs[:,iop],S[:,which-1])
    if norm(abs(dotvecs)-1)>eps:# We have at least one independent vector
        if which == 1: #Take the one that is most perpendicular to the earlier one
            found = True
            op2 = argmin(abs(dotvecs))
            newvec = rotvecs[:,op2]
        elif which == 2: # Third vector needs to be independent of first S vector, too. Use the first vector that is ind. of both
            for iop in range(parentlatt.nops):
                if isindependent(rotvecs[:,iop],S[:,0]) and isindependent(rotvecs[:,iop],S[:,1]):
                    found = True
                    newvec = rotvecs[:,iop]
                    break                                 
    return [found, newvec]
 
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
#        print; print iop; 
#        print 'symop R'; print parentlatt.symops[:,:,iop]
#        print 'symop m'; print dot(dot(inv(parentlatt.vecs[:,:]),parentlatt.symops[:,:,iop]),parentlatt.vecs[:,:])    
#        print 'det(m)', det( dot(dot(inv(parentlatt.vecs[:,:]),parentlatt.symops[:,:,iop]),parentlatt.vecs[:,:]) )  
#        print 'eigenvalues',evals
#        print 'eigenvectors'; print evecs #loop for printing only
        
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
            print 'Determinant zero'
            print lmat    
            return False
        mmat = dot(dot(inv(lmat),parentlatt.symops[:,:,iop]),lmat)
#        print 'mmat', iop
#        print mmat
        for i in range(3):
            for j in range(3):
                if abs(rint(mmat[i,j])-mmat[i,j])>1.0e-6:
#                    print iop, 'Symmetry failed for mmat[i,j]',mmat[i,j]
#                    print 'Cartesian operator' 
#                    print parentlatt.symops[:,:,iop] 
#                    print 'Cartesian Lattice'
#                    print lmat
#                    print 'transformed Cart operator' 
#                    print mmat                                           
                    return False #jumps out of subroutine
    return True #passes test

def nonDegen(vals):
     '''Tests whether a vector has one unique element.  If so, returns the index'''
     distinct = []
     if isreal(vals[0]) and not isequal(vals[0],vals[1]) and not isequal(vals[0],vals[2]):
         distinct.append(0)
     if isreal(vals[1]) and not isequal(vals[1],vals[0]) and not isequal(vals[1],vals[2]):
         distinct.append(1)
     if isreal(vals[2]) and not isequal(vals[2],vals[0]) and not isequal(vals[2],vals[1]):
         distinct.append(2)
     return distinct    

def matchDirection(vec,list):
    '''if vec parallel or antiparallel to any vector in the list, don't include it'''
    for vec2 in list:
        if isequal(abs(cosvecs(vec,vec2)),1.0):
            return True
    return False

def MT2mesh_one_n(MT,B):
        Q = dot(B.vecs,transpose(inv(MT)))
        N2 = rint((Nmesh/abs(det(MT)))**(1/3.0))
        ms = [N2,N2,N2]
#        print ms        
        Q[:,0] = Q[:,0]/ms[0]
        Q[:,1] = Q[:,1]/ms[1]
        Q[:,2] = Q[:,2]/ms[2]
        return Q

def MT2mesh_three_ns(MT,B):
        Q = dot(B.vecs,transpose(inv(MT)))
        mesh = svmesh(Nmesh/abs(det(MT)),Q)
        ms = mesh[0]
        print ms
        Q[:,0] = Q[:,0]/ms[0]
        Q[:,1] = Q[:,1]/ms[1]
        Q[:,2] = Q[:,2]/ms[2]
        return Q
     
##############################################################
########################## Script ############################

#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=200    

#print 'Target N kpoints', Nmesh

M = zeros((3,3),dtype = np_int)
S = zeros((3,3),dtype = np_float)
B = lattice()
A = lattice()
K = lattice()

##############BCT lattice
alat = 2*sqrt(5)
ca = 11/3.
clat = alat*ca
#B.vecs = matrix((  
#  [   -alat/2,  alat/2,   alat/2],
#  [   alat/2,  -alat/2,   alat/2],
#  [   clat/2,   clat/2,   -clat/2]
#  ), dtype=float)


#B.vecs = matrix((  #C axis along x !####
#  [   -clat/2,  clat/2,   clat/2],
#  [   alat/2,  -alat/2,   alat/2],
#  [   alat/2,   alat/2,   -alat/2]
#  ), dtype=float)

#print 'B vectors before inverse and transpose';print B.vecs
#B.vecs = trimSmall(transpose(inv(B.vecs)))
#############End BCT lattice

############## Any lattice

#crystal = [1,1,sqrt(2),90,90,120] # [a,b,c,alpha,beta,gamma]
crystal = [1,2,3,20,30,40] # [a,b,c,alpha,beta,gamma]
B.vecs = transpose(lattice_vecs(crystal))
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
A.Nmesh = Nmesh
print 'A vectors';print A.vecs
print 'Det of A', A.det
print 'Orth Defect of A', orthdef(A.vecs)
print 'Surf/vol of A', surfvol(A.vecs)

[A.symops,A.nops] = getGroup(A.vecs)
if A.nops != B.nops: 
    sys.exit('Number of operations different for A and B; stop')

print 'symmetry operations R of A\n'
testvecs = [];testindices = []
for k in range(A.nops):
    print; print k
    op = matrix(A.symops[:,:,k])
    print trimSmall(op)
    m = trimSmall(dot(dot(inv(A.vecs[:,:]), A.symops[:,:,k]),A.vecs[:,:])  ) 
    [vals,vecs]=eig(m); vecs = array(vecs)
#    print 'symop m'; print m
#
#    print 'det(m)', det(m)
    print 'eigen of m',vals
    print vecs
    print vecs[:,0]/abs(vecs[:,0])[nonzero(vecs[:,0])].min()
    print vecs[:,1]/abs(vecs[:,1])[nonzero(vecs[:,1])].min()
    print vecs[:,2]/abs(vecs[:,2])[nonzero(vecs[:,2])].min()
    
    #find operations with nondegenerate real eigenvalues
    print 'nonDegen', nonDegen(vals)
    for i in nonDegen(vals):
        if not matchDirection(transpose(vecs[:,i]),testvecs): #keep only unique directions    
            testvecs.append(vecs[:,i].real/abs(vecs[:,i])[nonzero(vecs[:,i])].min())
            testindices.append([k,i])
#print; print oplist;
print testvecs
print testindices
MT = zeros((3,3),dtype = np_int)


if len(testvecs) == 0:
    print 'No eigen directions'
#    MT = unconstrainedmin(B.vecs)
if len(testvecs) == 1:
    print 'Only 1 eigen direction'
    #Choose this one and the other two in the plane perpendicular to this. 
    #Since all symmetry operators will be diagonal in this mesh representaion
    #of eigenvectors of , 
    #
    MT[:,0] = testvecs[0]
    k = testindices[0][0]
    op = matrix(A.symops[:,:,k])
#    print trimSmall(op)
    m = trimSmall(dot(dot(inv(A.vecs[:,:]), A.symops[:,:,k]),A.vecs[:,:])  ) 
    [vals,vecs]=eig(m); vecs = array(vecs)
    print vecs
    otherindices = nonzero(array([0,1,2])-k)
    print otherindices[0]
    v1 = vecs[otherindices[0][0]]
    v2 = vecs[otherindices[0][1]]
    MT[:,1] = v1/abs(v1)[nonzero(v1)].min()
    #Make 3rd vector perp as possible to the other two 
    ur0 = dot(A.vecs,MT[:,0])/norm(dot(A.vecs,MT[:,0])) #unit vector in real space
    ur1 = dot(A.vecs,MT[:,1])/norm(dot(A.vecs,MT[:,1]))
    ur2 = cross(ur0,ur1)
    MT[:,2] = rint(dot(inv(A.vecs),ur2))
#    MT[:,2] = v2/abs(v2)[nonzero(v2)].min()  
    print 'MT from single operator';print MT
    Q2 = MT2mesh_three_ns(MT,B)
    if checksymmetry(Q2,B):
        SV = surfvol(Q2)
        print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'  
        K.vecs = Q2                
    else:
        print'Q fails symmetry'    

if len(testvecs) == 2:
    print 'Only 2 eigen directions'
    MT[:,0] = testvecs[0]
    MT[:,1] = testvecs[1]
    #Make 3rd vector perp as possible to the other two 
    ur0 = dot(A.vecs,MT[:,0])/norm(dot(A.vecs,MT[:,0])) #unit vector in real space
    ur1 = dot(A.vecs,MT[:,1])/norm(dot(A.vecs,MT[:,1]))
    ur2 = cross(ur0,ur1)
    MT[:,2] = rint(dot(inv(A.vecs),ur2))
    print 'MT from two eigen directions';print MT
    Q2 = MT2mesh_three_ns(MT,B)
    if checksymmetry(Q2,B):
        SV = surfvol(Q2)
        print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'  
        K.vecs = Q2                
    else:
        print'Q fails symmetry'  
                    
if len(testvecs) >= 3:
    testvecstrials = [list(x) for x in combinations(testvecs,3)]
    print testvecstrials    
    bestindex = -1 
    bestSV = 1000 
    for i,vecs in enumerate(testvecstrials):
        MT[:,0] = vecs[0]
        MT[:,1] = vecs[1]
        MT[:,2] = vecs[2]
        print 'MT from three eigen directions';print MT
    
        if not isequal(det(MT),0):
            Q2 = MT2mesh_one_n(MT,B)
            if checksymmetry(Q2,B):
                SV = surfvol(Q2)
                if SV<bestSV: bestSV = SV; bestindex = i; K.vecs = Q2
                print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'                  
            else:
                print'Q fails symmetry'
    print '___________ Best mesh ___________'
    print 'trial', bestindex
print K.vecs
K.det = abs(det(K.vecs))
print 'N of mesh', B.det/K.det
SV = surfvol(K.vecs)
print round(surfvol(K.vecs),4),round(orthdef(K.vecs),4),'SV of Q2,','OD' 


#sys.exit('stop')
#
#
#if len(testvecs)>0:
#    for vec in testvecs:
#        MT = zeros((3,3),dtype = np_int)
#        print;print '****************';print 'Cartesian direction to test as first superlattice vector '
#        print trimSmall(dot(A.vecs,vec))
#        print vec
#        MT[:,0] = rint(Nmesh**(1/3.0))*vec
##        print MT
#        S = array(dot(A.vecs,MT))
#        [found,newvec] = findNextVec(S,A,1)  #find next vector
#        if found:
#            print 'Found 2nd vector by symmetry'
#            S[:,1] =  transpose(newvec)
#            MT = dot(inv(A.vecs),S)
##            print MT
##            print S            
#            
#            [found,newvec] = findNextVec(S,A,2)
#            if found:
#                S[:,2] = transpose(newvec)
#                print 'Found 3rd vector by symmetry'
#                MT = dot(inv(A.vecs),S)
#                print S
##                print MT               
#            else: 
#                #minimize cost over ghi 
#                print 'Choosing third vector along an eigenvector, minimizing S/V'
#                MT = searchmin(S,A)
#                print MT    
#        else:
#            #minimize cost over defghi
#            print 'Choosing 2nd and 3rd vectors by minimizing S/V'
#            MT = searchmin(S,A)
#        S = dot(A.vecs,MT)
#        checksym = checksymmetry(S,A)
##        sys.exit('stop')
#        
#        
#        if checksym: 
#            K.vecs = transpose(inv(S))
#            checksymB = checksymmetry(K.vecs,B)
#            if checksymB: 
#    #                    print 'Obeys symmetry of lattice B:', checksymB               
#                print 'S';print S
#                print abs(det(S)/A.det),'Volume of superlattice'
#                print round(surfvol(S),2),round(orthdef(S),2),'SV of superlattice,','OD'  
#                print round(surfvol(K.vecs),2),round(orthdef(K.vecs),2),'SV of k-mesh,','OD'  
#                print 'M matrix abc;def;ghj';print trimSmall(transpose(dot(inv(A.vecs),S)))
#            else: 
#                print'Passed A symmetry, but failed B symmetry'
#        else:
#            print'Failed A symmetry'
#  
##        sys.exit('stop')
#else: #no symmetry directions
#    #minimize cost over defghi
#    print 'Finding all mesh vectors by minimizing S/V'
#    MT = searchmin(S,A)
#    S = dot(A.vecs,MT)
#    K.vecs = transpose(inv(S))
#    print 'S';print S
#    print abs(det(S)/A.det),'Volume of superlattice'
#    print round(surfvol(S),2),round(orthdef(S),2),'SV of superlattice','OD'  
#    print round(surfvol(K.vecs),2),round(orthdef(K.vecs),2),'SV of k-mesh','OD'  
#    print 'M matrix abc;def;ghj';print trimSmall(transpose(dot(inv(A.vecs),S)))   
#
##try random starting vector
##print '++++++++++++++++++++++++++'
##MT = searchminRandomStart(S,A)
##S = dot(A.vecs,MT)
##K.vecs = transpose(inv(S))
##print 'S';print S
##print abs(det(S)/A.det),'Volume of superlattice'
##print round(surfvol(S),2),round(orthdef(S),2),'SV of superlattice','OD'  
##print round(surfvol(K.vecs),2),round(orthdef(K.vecs),2),'SV of k-mesh','OD'  
##print 'M matrix abc;def;ghj';print trimSmall(transpose(dot(inv(A.vecs),S)))   
#
#
#           
#sys.exit('stop')
#
#
##if A.nops < 4: #has only 2 symm operations. 
#    
##else:
##    findmin(B,Nmesh)
##    
##    
##    minSV = 1000 #initialize
##    for a in range(maxint): # Diagonal elements of m can't be zero
##        for b in range(maxint):
##            for c in range(2*maxint):
##                if a==0 and b==0 and c==0:
##                    continue
##                #create first S vector
##                M[0,0]=a; M[0,1]=b; M[0,2]=c;
###                print; print [a,b,c],'a,b,c'
##                S[:,0] = dot(A.vecs,transpose(M[:,0]))
##                if norm(S[:,0])<eps:
##                    continue # all zeros, so go to next a,b,c
##                print '1st vector';print S[:,0]
##                #apply all symmetry operators, and find 2nd and 3rd vectors
##                S = fillS(S,A)
##                if abs(det(S))<eps:
##                    print 'S has zero det for ',[a,b,c]
##                    continue
##    #            print S
##                checksym = checksymmetry(S,A)
##                if checksym: 
##                    K.vecs = transpose(inv(S))
##                    checksymB = checksymmetry(K.vecs,B)
##                    if checksymB: 
##    #                    print 'Obeys symmetry of lattice B:', checksymB               
##                        print; print [a,b,c],'a,b,c'
##                        print 'S';print S
##                        print abs(det(S)/A.det),'Volume of superlattice'
##                        print round(surfvol(S),2),round(orthdef(S),2),'SV of superlattice','OD'  
##                        print round(surfvol(K.vecs),2),round(orthdef(K.vecs),2),'SV of k-mesh','OD'  
##                        print 'M matrix abc;def;ghj';print trimSmall(transpose(dot(inv(A.vecs),S)))
##                        if surfvol(K.vecs)<minSV:
##                            minSV = surfvol(K.vecs)
##                            bestK = K.vecs
##                            bestabc = [a,b,c]
##    print'------------------------------------------------'
##    print 'a,b,c',bestabc
##    print 'Best k mesh'; print bestK
##    print 'Best surface/volume:', round(minSV,2)
##
##print 'Done'