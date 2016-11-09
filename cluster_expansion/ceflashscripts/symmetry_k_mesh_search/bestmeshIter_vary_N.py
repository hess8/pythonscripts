#!/usr/bin/python
import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,  \
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
    #Fill a 1st brilloun zone with mesh points.  We will choose the 1st BZ to be that given by the parallepiped of (B0, B1, B2)
    #Since B = KM.  The first column of M determines the first column of B (B0) We run trial mesh points over a grid made by the maximum and minimum values of columns of M and the three directions 
    # of K.  The first row of M  gives the first k direction (first column of K)
    eps = 1e-4
    Kv = K.vecs
    Bv = B.vecs
    nBZpt = 0
    Binv = inv(Bv)
    print 'M in writekpts_vasp_M';print (M)
    print 'Kvecs in writekpts_vasp_M';print (Kv)
#    print 'transpose(Bvecs)in writekpts_vasp_M';print transpose(Bv)*100
    print 'det of M', det(M)    
    npts = -1
    ktryB = zeros((3,rint(det(M)*2)))# 
    kpts =  zeros((3,rint(det(M)*2)))
    #The rows of M determine how each vector (column) of K is used in the sum.    
    #The 1BZ parallelpiped must go from (0,0,0) to each of the other vertices 
    #the close vertices are at B1,B2,B3.  So each element of each row must be considered.
    #The far verictecs are at  for these three vectors taken in paris. 
    #To reach the diagonal point of the parallelpiped, 
    #which means that the sums of the rows must be part of the limits.
    #To reach the three far vertices (not the tip), we have to take the columns of M in pairs:, 
    #which means that we check the limits of the pairs among the elements of each row.
    #in other words, the limits on the search for each row i of (coefficients of grid basis vector Ki) are the partial sums
    #of the elements of each row:  min(0,a,b,c,a+b,a+c,b+c,a+b+c), max(0,a,b,c,a+b,a+c,b+c,a+b+c)
    Msums = zeros((3,8),dtype = int)

    for i in range(3):
        a = M[i,0]; b = M[i,1];c = M[i,2];
        Msums[i,0]=0; Msums[i,1]=a; Msums[i,2]=b;Msums[i,3]=c;
        Msums[i,4]=a+b; Msums[i,5]=a+c; Msums[i,6]=b+c;  Msums[i,7]=a+b+c
    ntry =0
    for i2 in range(amin(Msums[2,:])-1,amax(Msums[2,:])+1): #The rows of M determine how each vector (column) of M is used in the sum
        for i1 in range(amin(Msums[1,:])-1,amax(Msums[1,:])+1):
            for i0 in range(amin(Msums[0,:])-1,amax(Msums[0,:])+1):
                ntry += 1
                ktry = i0*Kv[:,0] + i1*Kv[:,1] + i2*Kv[:,2]              
                ktryB1 = trimSmall(dot(inv(Bv),transpose(ktry)))
               #test whether it is in 1st BZ.  Transform first to basis of B:
               #it's in the parallelpiped if its components are all less than one and positive             
                eps = 1e-4
                if min(ktryB1)>0-eps and max(ktryB1)<1-eps :
                    npts += 1
#                    print i0,i1,i2, trimSmall(ktryB1)
                    #translate to traditional 1BZ
                    for i in range(3):
                        if ktryB1[i]>0.5+eps: 
                            ktryB1[i] = ktryB1[i] - 1
                        if ktryB1[i]<-0.5+eps: 
                            ktryB1[i] = ktryB1[i] + 1
                    
                    #convert back to cartesian
                    ktry = trimSmall(dot(Bv,transpose(ktryB1)))
                    kpts[:,npts] = ktry
    npts = npts+1 #from starting at -1    
    print 'Grid points tested',ntry     
    print 'Points in 1BZ',npts
    if not areEqual(npts,rint(det(M))): 
        print det(M)
        sys.exit('Stop. Number of grid points in the 1BZ is not equal to det(M)')
    #Apply symmetry operations and see which are identical to others.  All in Cartesian coords
    kptssymm = zeros((3,npts))
    weights = zeros((npts),dtype = int)
    #record the first point
    kptssymm[:,0] = kpts[:,0]
    weights[0] = 1
    nksymm = 1
    
    for i in range(1,npts):
        kB = trimSmall(dot(inv(Bv),transpose(kpts[:,i])))
#        if areEqual(kB[0],0.5) or  areEqual(kB[1],0.5) or areEqual(kB[2],0.5)  :
#            print'Boundary point', kB
        #rotate
        found = False
        for iop in range(B.nops):
            krot = dot(B.symops[:,:,iop],kpts[:,i])
            kB2 = trimSmall(dot(inv(Bv),transpose(krot)))
#            if areEqual(kB[0],0.5) or  areEqual(kB[1],0.5) or areEqual(kB[2],0.5)  :
#                print kB2            
            #test whether it matches any we have saved. 
            for iksymm in range(nksymm):      
                if  areEqual(krot[0],kptssymm[0,iksymm]) and areEqual(krot[1],kptssymm[1,iksymm]) and areEqual(krot[2],kptssymm[2,iksymm]) :
#                    print 'Found equivalent point'
                    weights[iksymm] += 1
                    found = True # It better be equivalent to only one saved point
                    break
            if found: 
                break
        if not found:
            kptssymm[:,nksymm] = kpts[:,i]                
            weights[nksymm] += 1
            nksymm += 1  
#            print 'symm new point',nksymm  
    print 'Points in reduced 1BZ',nksymm 
    print 'Total weights',sum(weights)   
    print 'Vol BZ/ vol irredBZ', npts/float(nksymm)
    #convert to basis of B lattice vectors
    for i in range(nksymm):
        kptssymm[:,i] = trimSmall(dot(inv(Bv),transpose(kptssymm[:,i])))             
                                
#    #write POSCAR for vmd:  put B vectors in lattice, and kmesh in atomic positions
#    scale = 10       
#    poscar = open('POSCARk','w')
#    poscar.write('Cs I kpoints vs B'+'\n') #different sizes from this label
#    poscar.write('1.0\n')
#    for i in [0,1,2]:
#        poscar.write('%20.15f %20.15f %20.15f \n' % (scale*Bv[0,i], scale*Bv[1,i], scale*Bv[2,i])) 
#    poscar.write('1 %i\n' %npts)      
#    poscar.write('Cartesian\n')
#    poscar.write('0.0 0.0 0.0\n') 
#    for i in range(npts):
#        poscar.write('%20.15f %20.15f %20.15f \n' % (scale*kpts[0,i],scale*kpts[1,i],scale*kpts[2,i]))
#    poscar.close()
    
    #write POSCAR with irred BZ.  for vmd:  put B vectors in lattice, and kmesh in atomic positions
    scale = 10       
    poscar = open('POSCARkred','w')
    poscar.write('Cs I kpoints vs B'+'\n') #different sizes from this label
    poscar.write('1.0\n')
    for i in [0,1,2]:
        poscar.write('%20.15f %20.15f %20.15f \n' % (scale*Bv[0,i], scale*Bv[1,i], scale*Bv[2,i])) 
    poscar.write('1 %i\n' %nksymm)      
    poscar.write('Cartesian\n')
    poscar.write('0.0 0.0 0.0\n') 
    for i in range(nksymm):
        poscar.write('%20.15f %20.15f %20.15f %20.15f \n' % (scale*kptssymm[0,i],scale*kptssymm[1,i],scale*kptssymm[2,i], weights[i]))
    poscar.close()

    
    poscar = open('KPOINTS','w')
    poscar.write('BCH generated via bestmeshiter'+'\n') #different sizes from this label
    poscar.write('%i\n' % nksymm)
    poscar.write('Reciprocal lattice units\n')
    for i in range(nksymm):
        poscar.write('%20.15f %20.15f %20.15f      %i\n' % (kptssymm[0,i],kptssymm[1,i],kptssymm[2,i], weights[i]))
    poscar.close()
                    
#    sys.exit('stop')  
                
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

def bestmeshIter_vary_N(Blatt,Nmesh,path):
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
#    print'Symmetry operators in basis of B'
#    for i in range:
#        print B.msymops[:,:,i];print 
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
#    M0 = array([[2,   2,   2,],
#                    [2,   2,   -2],
#                    [-2,   2,   -2]])
#    
#      0.74050000    64.000 
#-5   1   -3   
#6   2   2   
#3   1   -3   
    M0 = array([[-5,   1, -3,],
                    [6,   2,   2],
                    [3,   1,   -3]])
    
#    for div in [256,128,64,32,16,8,4,2,1]:
#        print '\nDivisor',div
#        nMP = rint((Nmesh/div)**(1/3.0))
#        M = array([[nMP,0,0],[0,nMP,0],[0,0,nMP]]);
 

    for fac in [1,2,3,4,5,6,7,8]: 
        print '\nMultiplier',fac       
        
        M = fac*M0
#        M = array([[4,12,-4],
#                   [-11,4,-26],
#                   [-26,-4,-11]]);
        K = lattice();K.vecs = trimSmall(dot(B.vecs,inv(M)));K.det = abs(det(K.vecs)); K.Nmesh = B.det/K.det             
        print 'Number of points',det(M)
        print 'Check M'
        print M
        print 'Check K'
        print K.vecs 
        print 'Check B'
        print B.vecs
        print 'Check pf'
        print packingFraction(K.vecs) 
        #create a dir and prepare for vasp run
        newdir = str(K.Nmesh)
        newpath = path + newdir + '/'
        if not os.path.isdir(newpath):
            os.system('mkdir %s' % newpath)
        os.chdir(newpath)
        os.system ('cp %s* %s' % (vaspinputdir,newpath))
        os.system ('cp %sPOSCAR %s' % (path,newpath))  
        writekpts_vasp_M(newpath,B,M,K)
#             writekpts_vasp_pf(newpath,K.vecs,pf_maxpf,K.Nmesh)
        writejobfile(newpath)

#            print 'submitting job'            
        subprocess.call(['sbatch', 'vaspjob']) #!!!!!!! Submit jobs
        os.chdir(path)   
        

#'--------------------------------------------------------------------------------------------------------'
#'--------------------------------------------------------------------------------------------------------
