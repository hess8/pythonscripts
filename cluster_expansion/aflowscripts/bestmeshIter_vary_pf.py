import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms, intsymops, create_poscar, searchsphere

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil,real,unravel_index, outer, fmod, amin, amax

from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations
from pylab import frange

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

#def changewhichdual(Mv,B,run):
#    ''' delij: for i 1..9 and j 1..9, find the biggest negative cost change when varying 
#    i by -1,0,1 and j by -1,0,1.  If i = j, then we vary that single index by the first entry''' 
#    delij = zeros((9,9,2,2),dtype=float)#change in cost fo
#    oldcost = cost(Mv.reshape((3,3)),B,run)
#    for iv in range(9):
#        for jv in range(iv,9):
#            
#            if iv == jv:
#                for k,delm_k in enumerate([-1,1]):
#                    Mv[iv] += delm_k
#                    delij[iv,iv,k,k] = cost(Mv.reshape((3,3)),B,run) - oldcost
#                    Mv[iv] += - delm_k  
##                    print 'diag', iv,delij[iv,iv,k,k] 
#            else:
#                for k,delm_k in enumerate([-1,1]):
#                    for l,delm_l in enumerate([-1,1]):
#                        Mv[iv] += delm_k
#                        Mv[jv] += delm_l
#                        delij[iv,jv,k,l] = cost(Mv.reshape((3,3)),B,run) - oldcost
#                        Mv[iv] += -delm_k
#                        Mv[jv] += -delm_l 
##                        print iv,jv,delm_k,delm_l,delij[iv,jv,k,l]
#                
##    print 'min delij', min(delij)
#    print 'indices', unravel_index(argmin(delij), delij.shape)
#    print 'min delij'
#    print delij[unravel_index(argmin(delij), delij.shape)]
#
#    
##    print delij
#    return [unravel_index(argmin(delij), delij.shape),delij[unravel_index(argmin(delij), delij.shape)]]



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
        print 'Restart with different M'
        a = B.Nmesh**(1/3.0); c = 3;
        M = array([[-a+5, a/c , a/c],[a/c,-a,a/c],[a/c,a/c,-a-5]])
#        sys.exit('Stop')
    return trimSmall(M)
      


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
#            if run == 'minsvsym' and B.Nmesh/float(det(M))>1.2: M = M*2 # open up search space when when det(M) gets too low
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
#        Nscale =0*.05#.05; 
#        Ncost = Nscale * abs((B.det/det(kvecs))-B.Nmesh)/B.Nmesh 
#        shapescale = 0 * 0.01; shapecost = shapescale * surfvol(kvecs)
        cost = operr  #+ Ncost + shapecost

    return cost


def cost(M,B,run):
    if isequal(det(M),0):
        return 1000
    pftarget = B.pftarget
    
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
        cost = 1 * abs(pftarget - pf)/pftarget  + Ncost     
    elif run == 'minsvsym':
        Nscale =1*.05#.05; 
        Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        pfscale = 1 * 0.5; pfcost = pfscale * surfvol(K.vecs)
        symerr = symmetryError(K.vecs,B)
#        print symerr
#        cost = symerr *(1+Ncost)*(1+shapecost)
        cost = symerr  + Ncost + shapecost
    elif run == 'maxpfsym':
        Nscale =1*.2; #*.2;
        Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        pf = packingFraction(K.vecs)
        pfscale = 10 ; pfcost = pfscale * abs(pftarget - pf)/pftarget
        symerr = symmetryError(K.vecs,B)
#        print symerr          
#        cost = symerr *(1+Ncost)*(1+pfcost)
        cost = symerr  + Ncost + pfcost
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

def writekpts_vasp_pf(path,K,pf,Nmesh):
    '''Write mesh  to kpoints file, using integer division for cubic and fcc meshes'''   
    file1 = open(path +'KPOINTS','w')
    kpointsfile = []
    kpointsfile.append('%i kpoints for packing fraction pf=%6.4f\n' %(Nmesh,pf))
    kpointsfile.append('0 \n')   
    kpointsfile.append('Cartesian \n')
    for i in range(3):
        for j in range(3):
            kpointsfile.append('%18.12f' % K[j,i]) #transpose for Vasp input
        kpointsfile.append('\n')
#    kpointsfile.append('0.5 0.5 0.5\n' ) #shift
    kpointsfile.append('0.0 0.0 0.0\n' ) #shift
    file1.writelines(kpointsfile) 
    file1.close()
    return 

def writekpts_vasp_M(path,B,M,K):
    '''write out kpoints file with IBZKPTS format.  This will specify all the kpoints and their weights. 
    No shift is allowed for now'''
    #Fill a 1st brilloun zone with mesh points.  We will chose the 1st BZ to be that given by the parallepiped of (B0, B1, B2)
    #Since B = KM.  The first column of M determines the first column of B (B0) We run trial mesh points over a grid made by the maximum and minimum values of columns of M and the three directions 
    # of K.  The first row of M  gives the first k direction (first column of K).
    

    Kv = K.vecs
    Bv = B.vecs
    nBZpt = 0
    Binv = inv(Bv)
#    #Dummy set up Monkhorst Pack:
#    M = [[10,0,0],[0,10,0],[0,0,10]];
#    Kv = dot(B,inv(M))

#    
#    #end dummy
    print 'K vecsin writekpts_vasp_M';print (Kv)
#    print 'transpose(Bvecs)in writekpts_vasp_M';print transpose(Bv)*100
    print 'det of M', det(M)    
    npts = -1
    ktryB = zeros((3,int(det(M))))#!!!Change this *2 later when it's working
    kpts =  zeros((3,int(det(M))))
#    print size(ktryB)
#    for i2 in range(int(min(M[2,:])),int(max(M[2,:]))+10): #The rows of M determine how each vector (column) of M is used in the sum
#        for i1 in range(int(min(M[1,:])),int(max(M[1,:]))+10):
#            for i0 in range(int(min(M[0,:])),int(max(M[0,:]))+10):
#    print 'Min of M', amin(M)
#    print 'Max of M', amax(M)
#    imin = int(amin(M))
#    imax = int(amax(M))
#    for i2 in range(imin,imax): #The rows of M determine how each vector (column) of M is used in the sum
#        for i1 in range(imin,imax):
#            for i0 in range(imin,imax):
    [m0,m1,m2] = searchsphere(Kv)
    print [m0,m1,m2]
        
    for i2 in range(int(min(m2*M[2,:])),int(max(m2*M[2,:]))): #The rows of M determine how each vector (column) of M is used in the sum
        for i1 in range(int(min(m1*M[1,:])),int(max(m1*M[1,:]))):
            for i0 in range(int(min(m0*M[0,:])),int(max(m0*M[0,:]))):
                ktry = i0*Kv[:,0] + i1*Kv[:,1] + i2*Kv[:,2]
#                if i1 == 
#                print i0,i1,i2
#                print ktry
                
                ktryB1 = trimSmall(dot(inv(Bv),transpose(ktry)))
#                for i in range(npts):
#                    if ktryB1.all == ktryB[:,i].all:
#                        print 'Found duplicate'
#                        break
               #test whether it is in 1st BZ.  Transform first to basis of B:
               #it's in the parallelpiped if its components are all less than one and positive             
                eps = 1e-4
                if min(ktryB1)>=0+eps and max(ktryB1)<1+eps :
                    npts += 1
                    #translate to traditional 1BZ
                    for i in range(3):
                        if ktryB1[i]>0.5: 
                            ktryB1[i] = ktryB1[i] - 1
                        if ktryB1[i]<-0.5: 
                            ktryB1[i] = ktryB1[i] + 1
                    print i0,i1,i2, ktryB1
                    #convert back to cartesian
                    ktry = trimSmall(dot(Bv,transpose(ktryB1)))
                    kpts[:,npts] = ktry
                    
#                    ktryB[:,npts] = ktryB1
#                    kpts.append(ktryB1)
#                    print 'In 1BZ'
    npts = npts+1 #from starting at -1            
    print 'Points in 1BZ',npts
    #Apply symmetry operations and see which are identical to others.  All in Cartesian coords
    kptssymm = zeros((3,npts))
    weights = zeros((npts))
    #record the first point
    kptssymm[:,0] = kpts[:,0]
    weights[0] = 1
    nksymm = 1
    
    for i in range(1,npts):
        #rotate
        print i
        found = False
        for iop in range(B.nops):
            krot = dot(B.symops[:,:,iop],kpts[:,i])
            #test whether it matches any we have saved. 
            for iksymm in range(nksymm):      
                if  isequal(krot[0],kptssymm[0,iksymm]) and isequal(krot[1],kptssymm[1,iksymm]) and isequal(krot[2],kptssymm[2,iksymm]) :
                    print 'Found equivalent point'
                    weights[iksymm] += 1
                    found = True # It better be equivalent to only one saved point
                    break
   
            if found: 
                break
        if not found:
            kptssymm[:,nksymm] = kpts[:,i]                
            weights[nksymm] += 1
            nksymm += 1  
            print 'symm new point',nksymm  
    print 'Points in reduced 1BZ',nksymm 
    print 'Total weights',sum(weights)   
    print 'Vol BZ/ vol irredBZ', npts/float(nksymm)
                
                   
                                
    #write POSCAR for vmd:  put B vectors in lattice, and kmesh in atomic positions
    scale = 10       
    poscar = open('POSCARk','w')
    poscar.write('Cs I kpoints vs B'+'\n') #different sizes from this label
    poscar.write('1.0\n')
    for i in [0,1,2]:
        poscar.write('%20.15f %20.15f %20.15f \n' % (scale*Bv[0,i], scale*Bv[1,i], scale*Bv[2,i])) 
    poscar.write('1 %i\n' %npts)      
    poscar.write('Cartesian\n')
    poscar.write('0.0 0.0 0.0\n') 
    for i in range(npts):
        poscar.write('%20.15f %20.15f %20.15f \n' % (scale*kpts[0,i],scale*kpts[1,i],scale*kpts[2,i]))
    poscar.close()

                    
#    maxM = zeros((3,3),dtype = int)
#    minM = zeros((3,3),dtype = int)
#    N = zeros((3,3),dtype = int)
#    for i in range(3):
#        for j in range(3):
#            maxM[i,j] = int(max(M[i,j],0))
#            minM[i,j] = int(min(M[i,j],0))
#    print maxM
#    print minM
#    npts = 0
#    for i00 in range(minM[0,0],maxM[0,0]+1):
#        for i10 in range(minM[1,0],maxM[1,0]+1):
#            for i20 in range(minM[2,0],maxM[2,0]+1):
#                for i01 in range(minM[0,1],maxM[0,1]+1):
#                    for i11 in range(minM[1,1],maxM[1,1]+1):
#                        for i21 in range(minM[2,1],maxM[2,1]+1):
#                             for i02 in range(minM[0,2],maxM[0,2]+1):
#                                 for i12 in range(minM[1,2],maxM[1,2]+1):
#                                     for i22 in range(minM[2,2],maxM[2,2]+1):
#                                         N = [
#                                              [i00,i01,i02],
#                                              [i10,i11,i12],
#                                              [i20,i21,i22]
#                                              ]
#                                         pt = N
#                                         #convert to B basis
#                                         
#                                         


    sys.exit('stop')  
                
def writejobfile(path):
    '''read from a template in maindir, and put  (structure label) in job name'''
    file1 = open(path +'vaspjob','r')
    jobfile = file1.readlines()
    file1.close
    struct = path.split('/')[-3]
    pf = path.split('/')[-2]
    for i in range(len(jobfile)):
        jobfile[i]=jobfile[i].replace('myjob', struct+'_'+pf)
    file2 = open(path+'/'+'vaspjob','w')
    file2.writelines(jobfile) 
    file2.close()
    return

def bestmeshIter_vary_pf(Blatt,Nmesh,path):
    '''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer 3x3 matrix
    So K = B Inv(M).  Change M one element at a time to minimize the errors in symmetry and the cost in S/V and Nmesh '''
    
    ##############################################################
    ########################## Script ############################
    vaspinputdir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/vaspinput/'
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
#    B.pftarget = 0.7405 #default best packing fraction

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
    
    meshesfile = open('meshesfile','a')
    meshesfile.write('N target %i\n' % B.Nmesh)
    meshesfile.write('Format: pf then Nmesh then kmesh\n\n')    
    
    pflist = []
    for pftry in frange(pfB/2,0.75,0.005):
        print '\nPacking fraction target',pftry
        B.pftarget = pftry  
        pf_orth=0; pf_orth2fcc=0; sym_orth = False; sym_orth2fcc = False
#'--------------------------------------------------------------------------------------------------------'
#'--------------------------------------------------------------------------------------------------------'
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
#        print'Trial costs',ctest                           
        cbest = ctrials[argmin(ctest)]
#        print'Best c', cbest       
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
            M = rint(M)
            itersym = 0            
            symm = False
            while not symm and itersym <itermaxsym: 
                itersym += 1
                print 'Symmetry iteration', itersym, '-------'         
#                print 'Nmesh', abs(det(M)), 'packing', packingFraction(dot(B.vecs,inv(M)))
                M = rint(minkM(M,B))#; print'Mink reduced M'; print M    
                for iop in range(B.nops):
                    M = rint(findmin_i(M,B,iop))#need this to clean up numerical precision problems
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
                delcost = cost(M,B,type) - oldcost
        #write to files

#        meshesfile.write('Packing fraction target %f\n' % pftry)
        if symm and pf_maxpf not in pflist:
            pflist.append(pf_maxpf)
#            meshesfile.write('Packing fraction achieved %f\n' % pf_maxpf)
            meshesfile.write('%12.8f  %8.3f \n' % (pf_maxpf,K.Nmesh)) 
#            meshesfile.write('M\n')
#            for i in range(3):
#                for j in range(3):
#                    meshesfile.write('%i6' %M[i,j])
#                meshesfile.write('\n')
##            meshesfile.write('\n') 
                    
#            meshesfile.write('k mesh\n')
            for i in range(3):
                for j in range(3):
                    meshesfile.write('%18.12f' % K.vecs[i,j])
                meshesfile.write('\n')
            meshesfile.write('\n') 
            #create a dir and prepare for vasp run
            newdir = str(round(pf_maxpf,4))
            newpath = path + newdir + '/'
            if not os.path.isdir(newpath):
                os.system('mkdir %s' % newpath)
            os.chdir(newpath)
            os.system ('cp %s* %s' % (vaspinputdir,newpath))
            os.system ('cp %sPOSCAR %s' % (path,newpath))  
            writekpts_vasp_M(newpath,B,M,K)
            writekpts_vasp_pf(newpath,K.vecs,pf_maxpf,K.Nmesh)
            writejobfile(newpath)
            print 'Check M'
            print dot(inv(K.vecs),B.vecs) 
            print 'Check K'
            print K.vecs 
            print 'Check B'
            print B.vecs
            print 'Check pf'
            print packingFraction(K.vecs)  
#            print 'submitting job'            
            subprocess.call(['sbatch', 'vaspjob']) #!!!!!!! Submit jobs
            os.chdir(path)                      
        else:
            'do nothing'
#            meshesfile.write('Failed symmetry\n\n')     
    meshesfile.close()        
    
 #  Summary     
    pfs = [pfB]
    pftypes = ['B_latt']  
    ks  = [B.vecs/a] #one solutions is to simply divide B by an integer
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
