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
from numpy import ceil,real,unravel_index, outer, fmod

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

def changewhichdual(Mv,B,run):
    ''' delij: for i 1..9 and j 1..9, find the biggest negative cost change when varying 
    i by -1,0,1 and j by -1,0,1.  If i = j, then we vary that single index by the first entry''' 
    delij = zeros((9,9,2,2),dtype=float)#change in cost fo
    oldcost = cost(Mv.reshape((3,3)),B,run)
    for iv in range(9):
        for jv in range(iv,9):
            
            if iv == jv:
                for k,delm_k in enumerate([-1,1]):
                    Mv[iv] += delm_k
                    delij[iv,iv,k,k] = cost(Mv.reshape((3,3)),B,run) - oldcost
                    Mv[iv] += - delm_k  
#                    print 'diag', iv,delij[iv,iv,k,k] 
            else:
                for k,delm_k in enumerate([-1,1]):
                    for l,delm_l in enumerate([-1,1]):
                        Mv[iv] += delm_k
                        Mv[jv] += delm_l
                        delij[iv,jv,k,l] = cost(Mv.reshape((3,3)),B,run) - oldcost
                        Mv[iv] += -delm_k
                        Mv[jv] += -delm_l 
#                        print iv,jv,delm_k,delm_l,delij[iv,jv,k,l]
                
#    print 'min delij', min(delij)
    print 'indices', unravel_index(argmin(delij), delij.shape)
    print 'min delij'
    print delij[unravel_index(argmin(delij), delij.shape)]

    
#    print delij
    return [unravel_index(argmin(delij), delij.shape),delij[unravel_index(argmin(delij), delij.shape)]]


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
#            print 'newcost',newcost
            break
        else:
            M[bestindex[0],bestindex[1]] += bestindex[2]
#            print 'M altered in findmin_i';print M
            newcost = costi(M,B,iop)
#            print; print M;print newcost
        istep += 1
    if istep < maxsteps:
        'go on'
#        print 'Found minimum after %i steps' % istep
#        print 'Best M'; print M
    else:
        print 'Ended without minimum after %i steps' % istep
#        sys.exit('Stop')
    return trimSmall(M)

def gradM(Mv,B,run):
    '''Calculates the gradient in cost vs changing Mv by incrementing and decrementing each element'''
    gradM = zeros((9),dtype = float)
    for i in range(9):
        Mv[i] += 1;
        cost1 = cost(Mv.reshape((3,3)),B,run)
        Mv[i] += -2;
        cost2 = cost(Mv.reshape((3,3)),B,run)
        Mv[i] += 1  
        gradM[i] = (cost1-cost2)/2
    print 'gradM'; print gradM.reshape((3,3))
    return gradM

#def BFGSi(M,B,run):
#    '''BFGS method, but altered (BCH) to run over integers'''
#    Mv = M.flatten()
#    stepscale = 0.1 #Go only a fraction of step toward min location estimated by gradient
#    hessB = identity(9) #approx to hessian
#    delcost = -1.0
#    grad1 = gradM(Mv,B,run)
#    oldcost = cost(Mv.reshape((3,3)),B,run)
#    print 'cost',oldcost
#    while delcost < 0:
#        p = -dot(inv(hessB),grad1)
#        alpha = cost(Mv.reshape((3,3)),B,run)/norm(grad1)
#        print 'alpha',alpha;print 'p'; print p
#        s = rint(alpha*p*stepscale)
#        print 'det M',det(Mv.reshape((3,3)))
#        Mv = Mv + s
#        print 'new Mv'; print Mv.reshape((3,3))
#        print 'det M',det(Mv.reshape((3,3)))
#        grad2 = gradM(Mv,B,run)
#        y = grad2 - grad1
#        print 'y',y
##        print dot(B,outer(s,s))
##        hessB = hessB + outer(y,y)/inner(y,s) - dot(dot(hessB,outer(s,s)),hessB)/inner(s,dot(hessB,s))
#        newcost = cost(Mv.reshape((3,3)),B,run)
#        print 'cost',newcost
#        delcost = newcost - oldcost
#        print 'delcost',delcost
#        oldcost = newcost; grad1 = grad2 
##        sys.exit('stop')
#    M = Mv.reshape((3,3))
#    return M

#def findmin(M,B,run):
#    '''Dual version:  Finds minimum cost for the lattice by varying the integers of m, in the element that gives steepest descent
#    The 'run' indicates the cost function to use'''
#    maxsteps = 1000
#    istep = 1
#    M = minkM(M,B);print'Mink reduced M'; print M
#    Mv = M.flatten()
#   
#    while istep<maxsteps:
#
##        sys.exit('stop')
#        bestindex = changewhichdual(Mv,B,run)
#        print bestindex[0][0],bestindex[0][1], bestindex[0][2], bestindex[0][3],bestindex[1]
#        if bestindex[1]>=0:#found minimum
##            newcost = cost(M,B,run)
#            break
#        else:
#            if bestindex[0][0]==bestindex[0][1]: #only change this one index
#                Mv[bestindex[0][0]] += -(-1)**bestindex[0][2]#if index is 0, get -1; if index is 1, get 1
#            else:
#                Mv[bestindex[0][0]] += -(-1)**bestindex[0][2] 
#                Mv[bestindex[0][1]] += -(-1)**bestindex[0][3]
#        M = trimSmall(Mv.reshape((3,3)))
#        print 'New M'; print M
#        print 'newcost', cost(M,B,run)
#        if run == 'minsvsym' and B.Nmesh/float(det(M))>1.2: M = M*2 # open up search space when when det(M) gets too low
##            print; print M;print newcost
#        istep += 1
#    if istep < maxsteps:
#        print 'Found minimum after %i steps' % istep
##        print 'Best M'; print M
#        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));K.det = abs(det(K.vecs)); K.Nmesh = B.det/K.det             
#    else:
#        print 'Ended without minimum after maximum %i steps' % istep
#        sys.exit('Stop')
#    return [trimSmall(M),K]


def findmin(M,B,run): #normal routine for varying a single element at a time. 
    '''Finds minimum cost for the lattice by varying the integers of m, in the element that gives steepest descent
    The 'run' indicates the cost function to use'''
    maxsteps = 10000
    istep = 1
    M = minkM(M,B)#;print'Mink reduced M'; print M
    
    while istep<maxsteps:

#        sys.exit('stop')
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
#        print 'Found minimum after %i steps' % istep
#        print 'Best M'; print M
        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));K.det = abs(det(K.vecs)); K.Nmesh = B.det/K.det             
    else:
        print 'Ended without minimum after maximum %i steps' % istep
        sys.exit('Stop')
    return [trimSmall(M),K]

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
        cost = 1 * (0.7405 - pf)/0.7405  + Ncost     
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
        shapescale = 10 ; shapecost = shapescale * (0.7405 - packingFraction(K.vecs))/0.7405   
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
    pf_minsv = 0; pf_sv2fcc = 0; pf_maxpf = 0; pf_pf2fcc = 0; #kvecs_pf2fcc = identity(3)
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
        M = rint(dot(inv(kvecs_orth),B.vecs)).astype(int)
        print; print 'Try orth2FCC substitution.',
        kmesh2 = zeros((3,3),dtype = float)
        scale = 2/4**(1/3.0)
        kmesh2[:,0] = kvecs_orth[:,1]/scale + kvecs_orth[:,2]/scale
        kmesh2[:,1] = kvecs_orth[:,2]/scale + kvecs_orth[:,0]/scale
        kmesh2[:,2] = kvecs_orth[:,0]/scale + kvecs_orth[:,1]/scale   
        sym_orth2fcc = checksymmetry(kmesh2,B)
        if sym_orth2fcc:
            pf = packingFraction(kmesh2)
            print; print 'Packing fraction', pf, 'vs original B', pfB  
            if pf>pf_orth:
                M = rint(dot(inv(kmesh2),B.vecs)).astype(int)
                print 'M';print M
            else: 
                '    Packing fraction too small'            
            kvecs_orth2fcc = kmesh2
            pf_orth2fcc = pf
        else:
            print' It fails symmetry test'
#    sys.exit('stop')
    else:
#'--------------------------------------------------------------------------------------------------------'
#'--------------------------------------------------------------------------------------------------------'
   

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
        ctrials = [3]
        for c in ctrials:        
            M = array([[-a+2, a/c , a/c],[a/c,-a,a/c],[a/c,a/c,-a-2]])  #bcc like best avg pf on 50: 0.66
#            M = array([[2, a/c , a/c],[a/c,0,a/c],[a/c,a/c,-2]]) #fcc like best avg pf on 50: 0.59
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
        delcost = -1; lowcost = 1000 #initialize
####       while iternpf<itermaxnpf and delNPFcost <0 and abs(delNPFcost)>0.1 :
        while iternpf<itermaxnpf and delcost < -0.1:
            oldcost = cost(M,B,type)
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
                if symm:
                    newcost = cost(M,B,type)
                    if newcost - lowcost < 0: 
                        lowcost = newcost;
                        print'New lowcost',newcost              
                        pf_maxpf = packingFraction(K.vecs)
                        sym_maxpf = True
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
                        kvecs_pf2fcc = kmesh2
                        pf_pf2fcc = packingFraction(kmesh2)                        
                        Mtemp = rint(dot(inv(kmesh2),B.vecs)).astype(int)
                        if cost(Mtemp,B,type) < lowcost: 
                            lowcost = cost(M,B,type);print'New lowcost',lowcost
                            M = Mtemp                        
                            print 'M';print Mtemp
                            print 'Packing fraction', pf_pf2fcc, 'vs original B', pfB  

                            print;
                        else:
                            print 'Packing fraction too small' 
                    else:
                        print 'Fails to improve mesh'    
                delcost = cost(M,B,type) - oldcost
        
    #random test around best M:
    rrange = 3
    nrand = 300
    kmesh = trimSmall(dot(B.vecs,inv(M)))
    pfbest = packingFraction(kmesh)
    M2 = deepcopy(M)
    M3 = M2
    if 0.7405 - pfbest > 0.1:
        print 'Check %i random variations on M in range %i, vs best packing %f' % (nrand, rrange,pfbest)        
        for irun in range(nrand):
            if fmod(irun,100) == 0: print irun
            for i in range(3):
                for j in range(3):
#                        M3[i,j] = M3[i,j] + randint(-rrange,rrange+1)
                    M3[i,j] = M2[i,j] + randint(-rrange,rrange+1)
#            print;print 'Before scaling'; print M3
            M3 = rint(M3 * (B.Nmesh/abs(det(M3)))**(1/3.0))
#            print 'After scaling'; print M3              
            kmesh = trimSmall(dot(B.vecs,inv(M3)))
            if packingFraction(kmesh)>pfbest: #and not isequal(det(kmesh),0):
                [M3,K] = findmin(M3,B,type)
#                    for iop in range(B.nops):
#                        M3 = findmin_i(M3,B,iop)
                kmesh = trimSmall(dot(B.vecs,inv(M3)))
                if cost(M3,B,type) < lowcost and checksymmetry(kmesh,B):
                    lowcost = cost(M3,B,type)
                    pfbest = packingFraction(kmesh)
                    print "Rand M then search obeys sym"; print M3
                    print 'pf', packingFraction(kmesh)
                    print 'cost', cost(M3,B,type)
                    Nmesh = B.det/abs(det(kmesh))
                    print 'Nmesh', Nmesh, 'vs target', B.Nmesh                     
        
       
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
