import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint
from itertools import combinations
#
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

def findmin(M,B,run):
    '''Finds minimum cost for the lattice by varying the integers of m, in the element that gives steepest descent
    The 'run' indicates the cost function to use'''
    maxsteps = 1000
    istep = 1
    while istep<maxsteps:
        bestindex = changewhich(M,B,run)
        if bestindex[2]==0:#found minimum
            newcost = cost(M,B,run)
            break
        else:
            M[bestindex[0],bestindex[1]] += bestindex[2]
            newcost = cost(M,B,run)
        istep += 1
    if istep < maxsteps:
        print 'Found minimum after %i steps' % istep
#        print 'Best M'; print M
        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));K.det = det(K.vecs); K.Nmesh = B.det/K.det             
    else:
        print 'Ended without minimum after maximum %i steps' % istep
        sys.exit('Stop')
    return [trimSmall(M),K]

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
        Nscale =1*.05; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        shapescale = 1 * 0.5; shapecost = shapescale * surfvol(K.vecs)
        symerr = symmetryError(K.vecs,B)
#        print symerr
#        cost = symerr *(1+Ncost)*(1+shapecost)
        cost = symerr  + Ncost + shapecost
    elif run == 'maxpfsym':
        Nscale =1*.2; Ncost = Nscale * abs((B.det/K.det)-B.Nmesh)/B.Nmesh 
        shapescale = 1 * 0.5; shapecost = shapescale * (1/packingFraction(K.vecs))     
        symerr = symmetryError(K.vecs,B)
#        print symerr          
#        cost = symerr *(1+Ncost)*(1+shapecost)
        cost = symerr  + Ncost + shapecost
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
    pf_minsv = 0
    pf_sv2fcc = 0
    pf_maxpf = 0
    
       
    B.vecs = Blatt/2/pi  #Don't use 2pi constants in reciprocal lattice here.
    #############End BCT lattice
    eps = 1.0e-6

    B.Nmesh = Nmesh
    print 'B vectors';print B.vecs #
    #print 'B transpose'; print transpose(B.vecs)
    B.det = det(B.vecs)
    print 'Det of B', B.det
    print 'Orth Defect of B', orthdef(B.vecs)
    print 'Surf/vol of B', surfvol(B.vecs)
    pfB = packingFraction(B.vecs)
    print 'Packing fraction of B:', pfB  
    [B.symops,B.nops] = getGroup(B.vecs)
    printops_eigs(B)
    lattype = latticeType(B.nops)
    print 'Lattice type:', lattype
    
#    print 'MINK REDUCTION:'
#    B.vecs = transpose(mink_reduce(transpose(B.vecs), 1e-4)) #fortran routines use vectors as rows
#    print B.vecs
#    if not checksymmetry(B.vecs,B): sys.exit('Mink reduced lattice fails symmetry')
#    B.det = det(B.vecs)
#    print 'Det of B', B.det
#    print 'Orth Defect of B', orthdef(B.vecs)
#    print 'Surf/vol of B', surfvol(B.vecs)
#    pfB = packingFraction(B.vecs)
#    print 'Packing fraction of B:', pfB 
#    [B.symops,B.nops] = getGroup(B.vecs)
#    print 'Number of symmetry operations', B.nops
#    lattype = latticeType(B.nops)
#    print 'Mink lattice type:', lattype

    M = zeros((3,3),dtype=int)
    a = rint(Nmesh**(1/3.0)); f = int(Nmesh/a/a)
    #Check fcc-like mesh compatibility

    
#    if fcctype(B):

#    print 'fcc-like starting mesh'    
#    type = 'fcc' 
#    scale = 2*rint((Nmesh/4)**(1/3.0))
#    kmesh2 = zeros((3,3),dtype = float) 
#    kmesh2[:,0] = B.vecs[:,1]/scale + B.vecs[:,2]/scale
#    kmesh2[:,1] = B.vecs[:,2]/scale + B.vecs[:,0]/scale
#    kmesh2[:,2] = B.vecs[:,0]/scale + B.vecs[:,1]/scale 
##        print kmesh2
#    M = dot(inv(kmesh2),B.vecs)
#    print M    

    print;
    print 'Packing fraction minimization, starting with diagonal M'
    type = 'maxpf'  #revise these types  
    M[0,0]= a
    M[1,1]= a
    M[2,2]= f       
    [M,K] = findmin(M,B,type)
    print 'M ignoring symmetry:'; print M
    print 'Packing fraction:', packingFraction(K.vecs)
    print 'Nmesh', K.Nmesh, 'vs target', B.Nmesh  
    print 'Nearby mesh with required symmetry:'
    [M,K] = findmin(M,B,type+'sym')    
    print M
    kvecs_maxpf = K.vecs
    pf = packingFraction(kvecs_maxpf)
    sym_maxpf = checksymmetry(kvecs_maxpf,B)
    if sym_maxpf:
        print 'Packing fraction (maxpf):', pf, 'vs original B', pfB  
        print 'Nmesh', K.Nmesh, 'vs target', B.Nmesh        
        pf_maxpf = pf
    else: 
        print 'maxpf sym fail'; #status += 'maxpf sym fail;'; 
#    else:

    print;
    print 'S/V minimization, starting with diagonal M'
    type = 'minsv'  
    M[0,0]= a
    M[1,1]= a
    M[2,2]= f        
#    print M          
    [M,K] = findmin(M,B,type) 
    print 'M ignoring symmetry:'; print M
    print 'Packing fraction:', packingFraction(K.vecs)
    print 'Nmesh', K.Nmesh, 'vs target', B.Nmesh      
    print 'Nearby mesh with required symmetry:'          
    [M,K] = findmin(M,B,type+'sym')    
    print M
    kvecs_minsv = K.vecs
    pf = packingFraction(kvecs_minsv)   
    sym_minsv = checksymmetry(kvecs_minsv,B)
    sym_sv2fcc = False #initialize
    if sym_minsv:
        pf_minsv = pf
        print 'Packing fraction (minsv)', pf, 'vs original B', pfB  
        print 'Nmesh', K.Nmesh, 'vs target', B.Nmesh 
        print; print 'Try FCC-like substitution.',
        kmesh2 = zeros((3,3),dtype = float)
        scale = 2/4**(1/3.0)
        kmesh2[:,0] = K.vecs[:,1]/scale + K.vecs[:,2]/scale
        kmesh2[:,1] = K.vecs[:,2]/scale + K.vecs[:,0]/scale
        kmesh2[:,2] = K.vecs[:,0]/scale + K.vecs[:,1]/scale   
        sym_sv2fcc = checksymmetry(kmesh2,B)
        if sym_sv2fcc:
            print; print kmesh2
            pf = packingFraction(kmesh2)
            print 'Packing fraction', pf, 'vs original B', pfB  
            kvecs_sv2fcc = kmesh2
            pf_sv2fcc = pf
        else:
            print' It fails symmetry test'
                    
    else:
        print 'minsv sym fail'; #status += 'minsv sym fail;'
#
#    pfK = packingFraction(K.vecs)
#    if pfK < pfB: #didn't find better mesh; simply do Monkhorst-Pack with equal integers
#        meshtype = 'MHP' ; #status += 'MHPrevert;'
#        K.vecs = B.vecs/a
#        pfK = packingFraction(K.vecs)

    
#    Nratio = B.Nmesh/K.Nmesh #check to see if Nmesh is far enough off to correct by integer multiplication
#    if Nratio > 1.5:
#        mult = ceil(Nratio**(1/3.0))
#        M = M * mult        
##        mult = Nratio**(1/3.0)  
##        M = rint(M * mult)
#        print "Multiply M by", mult
#        run = 1
#        [M,K] = findmin(M,B,run)
#        print M
#        print 'Packing fraction:', round(packingFraction(K.vecs),4)  
        
#
#    print 'Final symmetry check:',checksymmetry(K.vecs,B)
    if not (sym_minsv or sym_sv2fcc or sym_maxpf):
        meshtype = 'B_latt_revert' ; #status += 'MHPrevert;'
        K.vecs = B.vecs/a; K.det = det(K.vecs); K.Nmesh = B.det/K.det
        pfmax = packingFraction(K.vecs)
    else:
        pfs = [pfB]
        pftypes = ['B_latt']
        ks  = [B.vecs/a]
        if sym_minsv:
            pfs.append(pf_minsv)
            pftypes.append('minsv')
            ks.append(kvecs_minsv)
        if sym_sv2fcc:
            pfs.append(pf_sv2fcc)
            pftypes.append('sv2fcc')
            ks.append(kvecs_sv2fcc)            
        if sym_maxpf:
            pfs.append(pf_maxpf)
            pftypes.append('maxpf')
            ks.append(kvecs_maxpf)             
        pfmax = max(pfs)
        meshtype = pftypes[argmax(pfs)]
        K.vecs = ks[argmax(pfs)]; K.det = det(K.vecs); K.Nmesh = B.det/K.det

#        print
#        print 'Final K mesh'; print K.vecs
#        print 'Final M'; print M
#        K.Nmesh = B.det/K.det
#        print 'N of mesh', B.det/K.det, 'vs target', B.Nmesh
#        print 'Packing fraction', pfK, 'vs original B', pfB            
#    else:
#        print'K mesh fails symmetry'        
#        sys.exit('Stop')
    
    return [K.vecs, K.Nmesh, B.Nmesh, lattype, pfB, pf_maxpf, pf_minsv, pf_sv2fcc, pfmax, meshtype, fcctype(B),status]
