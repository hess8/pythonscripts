import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms, intsymops

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil,real,unravel_index

from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations

def eigenvecfind(S,parentlatt):
    '''Find eigenvectors of symmetry operations (in cartesian basis).  A single eigenvector of eigenvalue 1 
    signifies a rotation operator. This returns a cartesian vector that corresponds to a rotation axis that is not in 
    the plane of the first two vectors of S'''
#    print 'S after array convert'; print S  
    eps = 1.0e-6
    
    for iop in range(parentlatt.nops): #loop for printing only
        evals,evecs = eig(parentlatt.symops[:,:,iop])            
    for iop in range(parentlatt.nops):
        evals,evecs = eig(parentlatt.symops[:,:,iop])      
        nOne = 0
        for i, eval in enumerate(evals):
            if abs(eval - 1) < eps:
                nOne += 1
                indexOne = i
        if nOne == 1:
            axis =  evecs[indexOne]
            unitaxis = axis/norm(axis) #unit vector
            if abs(cosvecs(cross(S[:,0],S[:,1]),axis))>eps: 
                #print'found independent third vector'
                return unitaxis
#    #print 'Found no 3rd vector through eigenvectors of symmetry operators'
#    sys.exit('stop')

def minkM(M,B):
    '''Transform to M on same lattice, but minked''' 
    kvecs = dot(B.vecs,inv(M))
#    print det(kvecs); print kvecs
    if det(kvecs) < 0.1: #do reduction in inverse space so we don't run into small number problems
        A = transpose(inv(B.vecs))
        S = dot(A,transpose(M))
        S = transpose(mink_reduce(transpose(S), 1e-4))
        M = transpose(dot(inv(A),S))
    else:
        kvecs = transpose(mink_reduce(transpose(kvecs), 1e-4)) #fortran routines use vectors as rows
        M = rint(dot(inv(kvecs),B.vecs))
    return trimSmall(M)
#
def changewhich_i(M,B,iop):
    bestgrad = 0
    bestdel = zeros((3,3),dtype=int)
    Mold = M
    oldcost = costi(Mold,B,iop)
#    print 'oldcost',oldcost
    bestindex=[-1,-1,0]#initialize
    for i in range(3):
        for j in range(3):
            M[i,j] += 1;delInc = costi(M,B,iop)-oldcost; M[i,j] += -1
            if delInc < 0 and delInc < bestgrad: bestindex = [i,j,1];bestgrad = delInc
            M[i,j] += -1;delDec = costi(M,B,iop)-oldcost;M[i,j] += 1;
            if delDec < 0 and delDec < bestgrad: bestindex = [i,j,-1];bestgrad = delDec
#            print i,j, delInc, delDec
    return bestindex

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

def grad(M,B,run):
    '''Explores changing M by two indices'''
    M = deepcopy(M)
#    bestgrad = 0
#    bestdel = zeros((3,3),dtype=int)
    delta = zeros((3,3),dtype = float)
    deltai = zeros((3,3),dtype = int)
    Mold = M
    oldcost = cost(Mold,B,run)
#    print 'oldcost',oldcost
    for i in range(3):
        for j in range(3):
            for inc in [-1,1]: #check increment and decrement
                M[i,j] += inc; 
                dcost = cost(M,B,run)-oldcost
                if dcost < delta[i,j]: 
                    delta[i,j] =dcost; deltai[i,j]=inc;
                M[i,j] += -inc
#    print 'delta';print delta; print 'deltai'; print deltai
    bestindex = unravel_index(argmin(delta), delta.shape)
#    print 'bestindex', bestindex
    M[bestindex]  += deltai[bestindex]
    # With M changed by one element, test all the indices to see the second index to change (could be same as first)
    delta2 = zeros((3,3),dtype = float)
    deltai2 = zeros((3,3),dtype = int)
    for i in range(3):
        for j in range(3):
            for inc in [-1,1]: #check increment and decrement:
                M[i,j] += inc; 
                dcost = cost(M,B,run)-oldcost
                if dcost < delta2[i,j]: 
                    delta2[i,j] =dcost; deltai2[i,j]=inc;
                M[i,j] += -inc    
    nextbest = unravel_index(argmin(delta2), delta2.shape); 
#    print 'delta2'; print delta2;print 'deltai2'; print deltai2
#    print "nextbest", nextbest
    return [bestindex, nextbest, deltai, deltai2]

def findmin_i(M,B,iop):
    '''Finds minimum cost for the lattice by varying the integers of m, in the element that gives steepest descent
    The 'run' indicates the cost function to use'''
    maxsteps = 200
    istep = 1
#    print 'M in findmin'; print M
    while istep<maxsteps:
        bestindex = changewhich_i(M,B,iop)
#        print 'bestindex',bestindex
        if bestindex[2]==0:#found minimum
            newcost = costi(M,B,iop)
            print 'newcost',newcost
            break
        else:
            M[bestindex[0],bestindex[1]] += bestindex[2]
#            print 'M altered in findmin_i';print M
            newcost = costi(M,B,iop)
#            print; print M;print newcost
        istep += 1
    if istep < maxsteps:
        print 'Found minimum after %i steps' % istep
#        print 'Best M'; print M
    else:
        print 'Ended without minimum after %i steps' % istep
#        sys.exit('Stop')
    return trimSmall(M)


def findmin(M,B,run):
    '''Finds minimum cost for the lattice by varying the integers of m, in the element that gives steepest descent
    The 'run' indicates the cost function to use'''
    maxsteps = 1000
    istep = 1
    M = minkM(M,B)#;print'Mink reduced M'; print M
    while istep<maxsteps:
        bestindex = changewhich(M,B,run)
        if bestindex[2]==0:#found minimum
#            newcost = cost(M,B,run)
            break
        else:
            M[bestindex[0],bestindex[1]] += bestindex[2]
#            newcost = cost(M,B,run)
            if run == 'minsvsym' and B.Nmesh/float(det(M))>1.2: M = M*2 # open up search space when when det(M) gets too low
#            print; print M;print newcost

        istep += 1
    if istep < maxsteps:
        print 'Found minimum after %i steps' % istep
#        print 'Best M'; print M
        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));K.det = abs(det(K.vecs)); K.Nmesh = B.det/K.det             
    else:
        print 'Ended without minimum after maximum %i steps' % istep
        sys.exit('Stop')
    return [trimSmall(M),K]

def gradM(M,B,run):
    '''Calculates the gradient in cost vs changing M by incrementing and decrementing each element'''
    gradM = zeros((3,3),dtype = float)
    oldcost = cost(M,B,run)
    bestindex=[-1,-1,0]#initialize
    for i in range(3):
        for j in range(3):
            M[i,j] += 1;
            cost1 = cost(M,B,run)
            M[i,j] += -2;
            cost2 = cost(M,B,run)
            M[i,j] += 1 #return to original state
            gradM[i,j] = (cost1-cost2)/2
    return gradM

#def findmin(M,B,run): #see above for other version
#    '''Finds minimum cost for the lattice by varying the integers of m, in the two elements that gives steepest descent
#    The 'run' indicates the cost function to use'''
#    maxsteps = 1000
#    istep = 1
#    while istep<maxsteps:
##        bestindex = changewhich(M,B,run)
#        [bestindex, nextbest, deltai, deltai2] = grad(M,B,run)
#        if deltai[bestindex]==0:#found minimum
#            break
#        else:
#            M[bestindex] += deltai[bestindex]
##            if random()>0.5: #don't do this every time to avoid cycles
#            M[nextbest] += deltai2[nextbest]
##            print 'bestindex, nextbest, deltai, deltai2', bestindex, nextbest; print deltai; print deltai2
##            print M
#            if run == 'minsvsym' and B.Nmesh/float(det(M))>1.2: M = M*2 # open up search space when when det(M) gets too low
##            print; print M;print newcost
#        istep += 1
#    if istep < maxsteps:
#        print 'Found minimum after %i steps' % istep
##        print 'Best M'; print M
#        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));K.det = abs(det(K.vecs)); K.Nmesh = B.det/K.det             
#    else:
#        print 'Ended without finding minimum, after the maximum %i steps' % istep
#        sys.exit('Stop')
#    return [trimSmall(M),K]

#def costsymflat(M,B):
#    M1 = deepcopy(M)
#    M1.shape = (3,3)
#    Kvecs = dot(B.vecs,inv(M1)) 
#    symerr = symmetryError(Kvecs,B)
#    return M1.flat   

#class cost2(object):  
#    def __init__(self):         
#        self.latt = []
#        self.symopi = []
#    def symflat(self,M):
#        M1 = deepcopy(M)
#        M1.shape = (3,3)
#        Kvecs = dot(self.latt.vecs,inv(M1)) 
#        symerr = symmetryError(Kvecs,self.latt)
#        print M1
#        print symerr
#        return symerr

def costi(M,B,iop):
    '''Here the symmetry cost is that due to only one symm operation,iop'''
    if det(M)<1: return 1000
    kvecs = dot(B.vecs,inv(M))
#    print 'iop in costi'; print B.symops[:,:,iop]
    mmat = trimSmall(dot(dot(inv(kvecs),B.symops[:,:,iop]),kvecs))
#    print 'mmat in costi';print mmat
    operr = 0.0
    for i in range(3):
        for j in range(3):
            if abs(rint(mmat[i,j])-mmat[i,j])>1.0e-4:
                operr += abs(rint(mmat[i,j])-mmat[i,j])
#                    print iop, 'Symmetry failed for mmat[i,j]',mmat[i,j]
#                    print 'Cartesian operator' 
#                    print parentlatt.symops[:,:,iop] 
#                    print 'Cartesian Lattice'
#                    print lmat
        Nscale =0*.05#.05; 
        Ncost = Nscale * abs((B.det/det(kvecs))-B.Nmesh)/B.Nmesh 
        shapescale = 0 * 0.01; shapecost = shapescale * surfvol(kvecs)
        cost = operr  + Ncost + shapecost

    return cost


def cost(M,B,run):
#    print M
    if isequal(det(M),0):
        return 1000
#        print 'Failed M'; print M
#        sys.exit('Det(M) = 0 in cost function. Stop')
#    else: 
#        K = lattice()
#        K.vecs = dot(B.vecs,inv(M));K.det = abs(det(K.vecs))
    K = lattice()
    K.vecs = dot(B.vecs,inv(M));K.det = abs(det(K.vecs))    
    if run == 'minsv':
        Nscale =1*1.0; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
#        cost = surfvol(K.vecs)*(1+Ncost)
        cost = surfvol(K.vecs) + Ncost
    elif run == 'maxpf':
#        K = lattice()
#        K.vecs = dot(B.vecs,inv(M));K.det = abs(det(K.vecs))
        Nscale =1*.5; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        pf = packingFraction(K.vecs)
#        cost = (1/pf)*(1+Ncost) 
        cost = 1/pf + Ncost     
    elif run == 'minsvsym':
        Nscale =1*.05#.05; 
        Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        shapescale = 1 * 0.5; shapecost = shapescale * surfvol(K.vecs)
        symerr = symmetryError(K.vecs,B)
#        print symerr
#        cost = symerr *(1+Ncost)*(1+shapecost)
        cost = symerr  + Ncost + shapecost
    elif run == 'maxpfsym':
        Nscale =1*.2; #*.2;
        Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        shapescale = 1 * 0.5; shapecost = shapescale * (1/packingFraction(K.vecs))     
        symerr = symmetryError(K.vecs,B)
#        print symerr          
#        cost = symerr *(1+Ncost)*(1+shapecost)
        cost = symerr  + Ncost + shapecost
    elif run == 'sym':   
        symerr = symmetryError(K.vecs,B)
#        print symerr
        cost = symerr 
    elif run == 'sym_sv':   
        symerr = symmetryError(K.vecs,B)
        shapescale = 1 * 0.5; shapecost = shapescale * surfvol(K.vecs)
        symerr = symmetryError(K.vecs,B)
        cost = symerr + shapecost        
    else:
       sys.exit('Cost type not found in cost function. Stop')      
    return(cost)  

def fcctype(B):
    '''See if an fcc-like mesh has the lattice symmetry'''
    latt2 = zeros((3,3),dtype = float)
    latt2[:,0] = B.vecs[:,1]/2 + B.vecs[:,2]/2
    latt2[:,1] = B.vecs[:,2]/2 + B.vecs[:,0]/2
    latt2[:,2] = B.vecs[:,0]/2 + B.vecs[:,1]/2   
    return checksymmetry(latt2,B)

def printops_eigs(B):
    print 'Number of symmetry operations', B.nops
    for j in range(B.nops):
        op = array(B.symops[:,:,j])
        [vals,vecs]=eig(op); vecs = array(vecs)
        print 'symop', j, 'egenvalues', vals
        print op
        print 'eigenvectors'; print vecs
        
def three_perp_eigs(A):
    testvecs = []; testindices = []
    svecs = zeros((3,3),dtype = float)
    for k in range(A.nops):
#        print; print k
        op = array(A.symops[:,:,k])
        [vals,vecs]=eig(op); vecs = array(vecs)      
        #find operations with nondegenerate real eigenvalues
#        print k, nonDegen(vals)
        for i in nonDegen(vals):
            if not matchDirection(transpose(vecs[:,i]),testvecs): #keep only unique directions    
                testvecs.append(real(vecs[:,i]))
                testindices.append([k,i])
    if len(testvecs) >= 3:
        testvecstrials = [list(x) for x in combinations(testvecs,3)]
#        print testvecstrials
        for trial in testvecstrials:
            #find superlattice
            for i in range(3):
#                print; print 'trial u',trial[i]
                svecs[i,:] = lattvec_u(A.vecs,trial[i])
#                print svecs[i,:]
#                print 'lattice m for lattice vector', dot(inv(A.vecs),transpose(svecs[i,:]))
#            print 'cosvecs', cosvecs(svecs[0,:],svecs[1,:]) , cosvecs(svecs[1,:],svecs[2,:]) , cosvecs(svecs[2,:],svecs[0,:])
#            print arenormal(svecs[0,:],svecs[1,:]) , arenormal(svecs[1,:],svecs[2,:]) , arenormal(svecs[2,:],svecs[0,:])           
            if arenormal(svecs[0,:],svecs[1,:]) and arenormal(svecs[1,:],svecs[2,:]) and arenormal(svecs[2,:],svecs[0,:]):
                S = transpose(array([svecs[0,:],svecs[1,:],svecs[2,:]]))
#                print 'found 3 perp'; print unique_anorms(S)
#                print S; print norm(S[:,0]); print norm(S[:,1]); print norm(S[:,2])
#                print unique_anorms(S).count(True); print A.lattype
                if A.lattype == 'Cubic' and unique_anorms(S).count(True) == 0:
                    return trimSmall(S) 
                if A.lattype == 'Tetragonal' and unique_anorms(S).count(True) <= 1:
                    return trimSmall(S) 
                if A.lattype == 'Orthorhombic': #and unique_anorms(S).count(True) == 3:
                    return trimSmall(S)
    sys.exit('error in three_perp_eigs')
        
def orthsuper(B):
    '''For lattice with orthogonal nonprimitive lattice vectors (cubic, tetragonal, orthorhombic),
    finds the simple orthorhombic superlattice with minimum s/v.'''
    # Find a set of three shortest lattice vectors that are perpendicular
    A = lattice()
    A.vecs = trimSmall(inv(transpose(B.vecs)))
#    print 'A'; print A.vecs
#    print 'det a', det(A.vecs)
    [A.symops,A.nops] = getGroup(A.vecs)
    A.lattype = latticeType(A.nops)
    
#    printops_eigs(A)
#    print 'transp A'; print transpose(A)    
    S = zeros((3,3),dtype = float)
    M = zeros((3,3),dtype = int)
    K = zeros((3,3),dtype = float)
#    S_orth = three_perp(A.vecs,B.lattype)
    print; S_orth =  three_perp_eigs(A)  
#    sys.exit('stop')  
    M = rint(transpose(dot(inv(A.vecs),S_orth))).astype(int)
    print 'M by finding 3 shortest perpendicular vectors';print M
    print 'det M', det(M)

    #starting mesh Q with 3 free directions. 
    Q = dot(B.vecs,inv(M))
#    print dot(Q[:,0],Q[:,1]),dot(Q[:,1],Q[:,2]),dot(Q[:,2],Q[:,0])
#    print norm(Q[:,0]),norm(Q[:,1]),norm(Q[:,2])
##        K = B.vecs/rint((B.Nmesh/det(Q))**(1/3.0))
#    if B.lattype == 'Tetrahedral':
#        #Find the free direction
#        q0 = norm(Q[:,0]); q1 = norm(Q[:,1]); q2 = norm(Q[:,2]); 
#        if isequal(q1,q2)): nmesh 
    print 'mesh numbers'; 
    [n0,n1,n2] = svmesh(B.Nmesh/abs(det(M)),Q)[0]
    print [n0,n1,n2]
    K[:,0] = Q[:,0]/n0; K[:,1] = Q[:,1]/n1; K[:,2] = Q[:,2]/n2
#    print K
    Nmesh = B.det/abs(det(K))
    if checksymmetry(K,B):
        pf = packingFraction(K)
        print 'Packing fraction (orthmesh)', pf, 'vs original B', packingFraction(B.vecs)  
        print 'Nmesh', Nmesh, 'vs target', B.Nmesh 
    else:
        sys.exit('Symmetry failed in orthsuper')
    return [K,pf,True]
    
def bestmeshIter(Blatt,Nmesh):
    '''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer 3x3 matrix
    So K = B Inv(M).  Change M one element at a time to minimize the errors in symmetry and the cost in S/V and Nmesh '''
    
    ##############################################################
    ########################## Script ############################
   
    M = zeros((3,3),dtype = int)
    S = zeros((3,3),dtype = fprec)
    B = lattice()
    A = lattice()
    K = lattice()
    status = ''
    pf_minsv = 0; pf_sv2fcc = 0; pf_maxpf = 0; pf_pf2fcc = 0
    sym_maxpf = False;  sym_sv2fcc = False; sym_minsv = False; sym_pf2fcc = False
    a = rint(Nmesh**(1/3.0)); f = int(Nmesh/a/a)
    print 'Target mesh number', Nmesh
       
    B.vecs = Blatt/2/pi  #Don't use 2pi constants in reciprocal lattice here
    #############End BCT lattice
    eps = 1.0e-6

    B.Nmesh = Nmesh
    print 'B vectors (differ by 2pi from traditional)';print B.vecs #
    #print 'B transpose'; print transpose(B.vecs)
    B.det = det(B.vecs)
    print 'Det of B', B.det
    print 'Orth Defect of B', orthdef(B.vecs)
    print 'Surf/vol of B', surfvol(B.vecs)
    pfB = packingFraction(B.vecs)
    print 'Packing fraction of B:', pfB  
    [B.symops,B.nops] = getGroup(B.vecs)
    B.msymops = intsymops(B) #integer sym operations in B basis
#    printops_eigs(B)
    B.lattype = latticeType(B.nops)
    print 'Lattice type:', B.lattype
    A = lattice()
    A.vecs = trimSmall(inv(transpose(B.vecs)))
    [A.symops,A.nops] = getGroup(A.vecs)    
    A.msymops = intsymops(A)
    print 'Real space lattice A'; print A.vecs
    print 'Det A', det(A.vecs)
    pfA = packingFraction(A.vecs)
    print 'Packing fraction of A:', pfA    
    
    pf_orth=0; pf_orth2fcc=0; sym_orth = False; sym_orth2fcc = False
    if B.lattype in ['Orthorhombic', 'Tetragonal','Cubic']:
        cbest = '' #need for passing other structures' results to main program 
        [kvecs_orth,pf_orth,sym_orth] = orthsuper(B)
        print; print 'Try orth2FCC substitution.',
        kmesh2 = zeros((3,3),dtype = float)
        scale = 2/4**(1/3.0)
        kmesh2[:,0] = kvecs_orth[:,1]/scale + kvecs_orth[:,2]/scale
        kmesh2[:,1] = kvecs_orth[:,2]/scale + kvecs_orth[:,0]/scale
        kmesh2[:,2] = kvecs_orth[:,0]/scale + kvecs_orth[:,1]/scale   
        sym_orth2fcc = checksymmetry(kmesh2,B)
        if sym_orth2fcc:
            M = rint(dot(inv(kmesh2),B.vecs)).astype(int)
            print; print 'M';print M
            pf = packingFraction(kmesh2)
            print; print 'Packing fraction', pf, 'vs original B', pfB  
            kvecs_orth2fcc = kmesh2
            pf_orth2fcc = pf
        else:
            print' It fails symmetry test'
#    sys.exit('stop')
    else:
#'--------------------------------------------------------------------------------------------------------'
#'--------------------------------------------------------------------------------------------------------'
#        print 'MINK REDUCTION:'
#        B.vecs = transpose(mink_reduce(transpose(B.vecs), 1e-4)) #fortran routines use vectors as rows
#        print B.vecs
#        if not checksymmetry(B.vecs,B): sys.exit('Mink reduced lattice fails symmetry')
#        B.det = det(B.vecs)
#        print 'Det of B', B.det
#        print 'Orth Defect of B', orthdef(B.vecs)
#        print 'Surf/vol of B', surfvol(B.vecs)
#        pfB = packingFraction(B.vecs)
#        print 'Packing fraction of B:', pfB 
#        [B.symops,B.nops] = getGroup(B.vecs)
#        print 'Number of symmetry operations', B.nops
#        B.lattype = latticeType(B.nops)
#        print 'Mink lattice type:', B.lattype
     

#        if B.lattype in ['Hexagonal','Monoclinic']:
#            M = zeros((3,3),dtype=int)
##             print 'Starting with diagonal M'
##             M[0,0]= a
##             M[1,1]= a
##             M[2,2]= f
#            c=3        
#            M = array([[-a, a/c , a/c],[a/c,-a,a/c],[a/c,a/c,-a]])
#     
#        else:
        M = zeros((3,3),dtype=int)
        ctest = []
        type = 'maxpfsym'; print type
        ctrials = [3,1000]
        for c in ctrials:        
            M = array([[-a+2, a/c , a/c],[a/c,-a,a/c],[a/c,a/c,-a-2]])  
            M = rint(M * (B.Nmesh/abs(det(M)))**(1/3.0))
            print 'Start mesh trial'; print M              
            [M,K] = findmin(M,B,type)
            print 'Test trial M'; print M

            ctest.append(cost(M,B,type))
        print'Trial costs',ctest                           
        cbest = ctrials[argmin(ctest)]
        print'Best c', cbest
            
            
        
     #
     #        print 'S/V minimization, start diag a,w/ random off diag'
     #        for i in range(3):
     #            for j in range(3):
     #                if i ==j:
     #                    M[i,j] = a
     #                else:
     #                    M[i,j] = 3 
     #                    M[i,j] = rint(a/3)
#        print M  
     #        ma = zeros((3,3,A.nops),dtype = float)
     #        ms = zeros((3,3,A.nops),dtype = float)
     #        for iop in range(A.nops): ma[:,:,iop] = trimSmall(dot(dot(inv(A.vecs),A.symops[:,:,iop]),A.vecs))
     #        ms[:,:,iop] = trimSmall(dot(dot(inv(S),ma[:,:,iop]),S)) #starting values
     #        type = 'minsv';print type
     
     #        type = 'maxpf'; print type
     #        [M,K] = findmin(M,B,type) 
     #        print 'M ignoring symmetry:'; print M
       
        iternpf = 0
        itermaxnpf = 10
        itermaxsym = 5
     #        bestpf = 100
#        NPFcost = 100;delNPFcost = -1 #initial values
        type = 'maxpfsym'; print type
        lowcost = 1000;
#        while iternpf<itermaxnpf and delNPFcost <0 and abs(delNPFcost)>0.1 :
        while iternpf<itermaxnpf and cost(M,B,type) <= lowcost:
        #            NPFcost = cost(M,B,'maxpf')
#            delNPFcost = (NPFcost-oldNPFcost)/NPFcost :
            '''Here we let M vary in the search, but record pf and kvecs when we find min cost'''
#            print 'cost(N,PF):', cost(M,B,type)
     #        while not symm and and iternpf<itermax:
            M = rint(M * (B.Nmesh/abs(det(M)))**(1/3.0))
            print 'Scaled M';print M
            iternpf += 1
            print 'Iteration',type,iternpf, '**********'
            [M,K] = findmin(M,B,type) 
            print M
            itersym = 0            
            symm = False
            while not symm and itersym <itermaxsym: 
                itersym += 1
                print 'Symmetry iteration', itersym, '-------'         
#                if abs((B.det/det(K.vecs))-B.Nmesh)/B.Nmesh > 0.15: #how far off from target N
#                    M = rint(M * (B.Nmesh/abs(det(M)))**(1/3.0))
#                    print 'Scaled M';print M
                print 'Nmesh', abs(det(M)), 'packing', packingFraction(dot(B.vecs,inv(M)))
                M = minkM(M,B)#; print'Mink reduced M'; print M    
                for iop in range(B.nops):
                    M = findmin_i(M,B,iop)
                    if abs(det(M)-B.Nmesh)/B.Nmesh > 0.15: #how far off from target N
                        M = rint(M * (B.Nmesh/abs(det(M)))**(1/3.0))
                        print 'Scaled M';print M                    
                K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));K.det = abs(det(K.vecs)); K.Nmesh = B.det/K.det                                       
                symm = checksymmetry(K.vecs,B)
                print 'Symmetry check', symm
                print M
                if symm:
                    if cost(M,B,type) < lowcost: 
                        lowcost = cost(M,B,type);
                        print'New lowcost',lowcost              
                        pf_maxpf = packingFraction(K.vecs)
                        sym_maxpf = symm
                        kvecs_maxpf = K.vecs
     #                    if pf_maxpf<bestpf: bestpf = pf_maxpf; bestM = M
                    print 'Packing fraction', pf_maxpf, 'vs original B', pfB  
                    print 'Nmesh', K.Nmesh, 'vs target', B.Nmesh 
                    print; print 'Try FCC-like substitution.'
                    kmesh2 = zeros((3,3),dtype = float)
                    scale = 2/4**(1/3.0)
                    kmesh2[:,0] = K.vecs[:,1]/scale + K.vecs[:,2]/scale
                    kmesh2[:,1] = K.vecs[:,2]/scale + K.vecs[:,0]/scale
                    kmesh2[:,2] = K.vecs[:,0]/scale + K.vecs[:,1]/scale 
        #            M = rint(dot(inv(kmesh2),B.vecs)).astype(int) #set this for maxpf run  
                    if checksymmetry(kmesh2,B):
                        sym_pf2fcc = True
                        Mtemp = rint(dot(inv(kmesh2),B.vecs)).astype(int)
                        if cost(Mtemp,B,type) < lowcost: 
                            lowcost = cost(M,B,type);print'New lowcost',lowcost  
                            M = Mtemp                   
                            print; print 'M';print M            
            #            print; print kmesh2
                        pf_pf2fcc = packingFraction(kmesh2)
     #                        if pf_sv2fcc<bestpf: bestpf = pf_sv2fcc; bestM = M
                        print 'Packing fraction', pf_pf2fcc, 'vs original B', pfB  
                        kvecs_pf2fcc = kmesh2
                    else:
                        print' It fails symmetry test'                     
#            NPFcost = cost(M,B,'maxpf')
#            delNPFcost = (NPFcost-oldNPFcost)/NPFcost
#            print 'Change in N,PF cost:', round(delNPFcost,4)
#            oldNPFcost = NPFcost
     #            print NPFcost

           
#        sys.exit('stop')       
 
       
 #  Summary     
    pfs = [pfB]
    pftypes = ['B_latt']  
    ks  = [B.vecs/a]   
    if not (sym_minsv or sym_sv2fcc or sym_maxpf or pf_pf2fcc or sym_orth or sym_orth2fcc):
         meshtype = 'B_latt_revert' ; #status += 'MHPrevert;'
         K.vecs = B.vecs/a; K.det = abs(det(K.vecs)); K.Nmesh = abs(B.det/K.det)
         pfmax = packingFraction(K.vecs)
    else:     

         if sym_orth:
             pfs.append(pf_orth)
             pftypes.append('orth')
             ks.append(kvecs_orth)
         if sym_orth2fcc:
            pfs.append(pf_orth2fcc)
            pftypes.append('orth2fcc')
            ks.append(kvecs_orth2fcc)          
#        if sym_minsv:
#            pfs.append(pf_minsv)
#            pftypes.append('minsv')
#            ks.append(kvecs_minsv)
#        if sym_sv2fcc:
#            pfs.append(pf_sv2fcc)
#            pftypes.append('sv2fcc')
#            ks.append(kvecs_sv2fcc)            
         if sym_maxpf:
             pfs.append(pf_maxpf)
             pftypes.append('maxpf')
             ks.append(kvecs_maxpf)    
         if sym_pf2fcc:
             pfs.append(pf_pf2fcc)
             pftypes.append('pf2fcc')
             ks.append(kvecs_pf2fcc)                    
    pfmax = max(pfs)
    meshtype = pftypes[argmax(pfs)]
    K.vecs = ks[argmax(pfs)]; K.det = abs(det(K.vecs)); K.Nmesh = B.det/K.det
#    return [K.vecs, K.Nmesh, B.Nmesh, B.lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_minsv, pf_sv2fcc, pfmax, meshtype, fcctype(B),status]

    return [K.vecs, K.Nmesh, B.Nmesh, B.lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype(B),cbest,status]
