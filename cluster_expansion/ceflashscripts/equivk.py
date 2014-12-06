#!/usr/bin/python
import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms, intsymops, create_poscar, searchsphere

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil,real,unravel_index, outer, fmod, amin, amax, sign

#from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations
from pylab import frange


from ceroutines import readfile,readstructs,readlatt,scaleposcar, getscale, getscalePOSCAR, getline, getL, \
    writekpts_cubic_n,writekpts_fcc_n,writefile, makestr2poscar

def writekpts_vasp_M(path,Bv,M,nops,symmops):
    '''write out kpoints file with IBZKPTS format.  This will specify all the kpoints and their weights. 
    No shift is allowed for now.  The symmops are those of B'''
    #Fill a 1st brilloun zone with mesh points.  We will choose the 1st BZ to be that given by the parallepiped of (B0, B1, B2)
    #Since B = KM.  The first column of M determines the first column of B (B0) We run trial mesh points over a grid made by the maximum and minimum values of columns of M and the three directions 
    # of K.  The first row of M  gives the first k direction (first column of K)
    eps = 1e-4
    M = M * sign(det(M)) # want positive determinants
    Kv = dot(B,inv(M))
    nBZpt = 0
    Binv = inv(Bv)
    print 'M in writekpts_vasp_M';print (M)
    print 'Kvecs in writekpts_vasp_M';print (Kv)
#    print 'transpose(Bvecs)in writekpts_vasp_M';print transpose(Bv)*100
    print 'det of M', det(M)    
    npts = -1
    ktryB = zeros((3,rint(det(M)*2))) 
    kpts =  zeros((3,rint(det(M)*2))) #all the kpoints in the entire 1BZ
    
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
    if not isequal(npts,rint(det(M))): 
        print det(M)
        sys.exit('Stop. Number of grid points in the 1BZ is not equal to det(M)')
    #Apply symmetry operations and see which are identical to others.  All in Cartesian coords
    kptssymm = zeros((3,npts)) #the kpoints in irreducible 1BZ
    weights = zeros((npts),dtype = int)
    #record the first point
    print '0'; print kpts[:,0]
    kptssymm[:,0] = kpts[:,0]
    weights[0] = 1
    nksymm = 1
    
    for i in range(1,npts): #test all 
        print; print i;print kpts[:,i]
        kB = trimSmall(dot(inv(Bv),transpose(kpts[:,i])))#now in the basis of B vectors
        print kB, 'in recip cords'
#        if isequal(kB[0],0.5) or  isequal(kB[1],0.5) or isequal(kB[2],0.5)  :
#            print'Boundary point', kB
        #rotate
        found = False
        for iop in range(nops):
            krot = dot(symops[:,:,iop],kpts[:,i])
            kB2 = trimSmall(dot(inv(Bv),transpose(krot)))
#            if isequal(kB[0],0.5) or  isequal(kB[1],0.5) or isequal(kB[2],0.5)  :
#                print kB2            
            #test whether it matches any we have saved. 
            for iksymm in range(nksymm):      
                if  isequal(krot[0],kptssymm[0,iksymm]) and isequal(krot[1],kptssymm[1,iksymm]) and isequal(krot[2],kptssymm[2,iksymm]) :
                    print 'Found equivalent point',iksymm;print kptssymm[:,iksymm]
                    weights[iksymm] += 1
                    found = True # It better be equivalent to only one point in irreducible 1BZ
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
    jobfile = readfile(path +'vaspjob')
#    print path.split('/'); sys.exit()
    label = path.split('/')[-2]
    for i in range(len(jobfile)):
        jobfile[i]=jobfile[i].replace('myjob',label)
    file2 = open(path+'/'+'vaspjob','w')
    file2.writelines(jobfile) 
    file2.close()
    return

#def writeincar(path):
#    poscar = open('KPOINTS','w')
#    poscar.write('BCH generated via bestmeshiter'+'\n') #different sizes from this label
#    poscar.write('%i\n' % nksymm)
#    poscar.write('Reciprocal lattice units\n')
#    for i in range(nksymm):
#        poscar.write('%20.15f %20.15f %20.15f      %i\n' % (kptssymm[0,i],kptssymm[1,i],kptssymm[2,i], weights[i]))
#    poscar.close()
#
#def getngxf(path):
# ==================================================================================
# ==================================================================================
# ==================================================================================
''' 
Here we read lattice in a 'scale free' manner, i.e. every cubic lattic is 100, 010, 001. 
The real lattice is A = PL, where P is the parent lattice and L is an integer matrix.  
We are forming a mesh K in the form B = KM, then scale it for use.   

IF the K mesh is cubic, then M is nB (again scale free), where n is a multiple (m) of
the volume factor det(L), and B = trans(inv(A))

If the Kmesh is fcc, then...



'''

#atomic = 'Al:Al'
#atomic = 'Si:Si'
#maindir = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6_encut500/'
#maindir = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6.prec.accurate/'

#maindir = '/fslhome/bch/cluster_expansion/sisi/equivk_accurate/'
#maindir = '/fslhome/bch/cluster_expansion/alal/cubic_al/equivk_c1-6_accurate/'
maindir = '/fslhome/bch/cluster_expansion/alal/equivk_f-16.tetra.noBlochl/'
#maindir = '/fslhome/bch/cluster_expansion/cucu/equivk_accurate/'
#maindir = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6.tetra/'
#maindir = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6.sigma.02/'
#enumfile = maindir + 'struct_enum.in.si'
#enumfile = maindir + 'struct_enum.in.cub' 
enumfile = maindir + 'struct_enum.in.fcc' 
#structfile = '../c1_6.dat'
#structfile = '../f1.dat'
structfile = '../f1_6.dat'
#finalDir = '/fslhome/bch/cluster_expansion/alir/enumtest/structs.myk/'
#finalDir = '/fslhome/bch/cluster_expansion/alir/enumtest/structs.cubicmesh/'
#finalDir = '/fslhome/bch/cluster_expansion/alir/enumtest/structs.cubictest/'
#finalDir = '/fslhome/bch/cluster_expansion/alir/enumtest/structs.fccmesh/'
finalDir = maindir + 'structs.cubmesh/'
#finalDir = maindir + 'structs.fccmesh/'
nmax = 45
print finalDir.split('.')
if 'cub' in finalDir.split('.')[-1]:
    meshtype = 'cub'
    multlist = []
#    multlist = [2,3,4,5,6,7,8]
    multlist = [2,3,4,5,6,8,10,12,13,14,16,18,19,20,21,22,23,24] #this is m
#    multlist = [2,3,4,5,6,8,10,12,13,14,16,18,19,20,21,22,23,24,26,28,32,36,38,40,42,44] #this is m
#    multlist = [40,42,44] #this is m
#    multlist = [2]    

elif 'fcc' in finalDir.split('.')[-1]:
    meshtype = 'fcc'
#    multlist = [2,3]
#    multlist = [2,3,4,5]  #5*4^(1/3) =~ 8, for fcc scaling
#    multlist = [2,3,4,5,6,8,10,12,13,14,16,18,19,20,21,22,23,24,26,28]  #44/4^(1/3) = 28
    multlist = [2,3,4,5,6,8,10,12]  #44/4^(1/3) = 28
else:
    sys.exit('finalDir must contain cub or fcc substring')
    
if not os.path.isdir(finalDir): os.system('mkdir %s' % finalDir)
vaspinputdir = maindir + 'vaspinput/'  
   
structchar = getline(0,enumfile).split('.')[-1][0]


fccform = array([[0,1,1],[1,0,1],[1,1,0]])
bccform = array([[1,1,-1],[1,-1,1],[-1,1,1]])
bccsq = dot(bccform,bccform) #

#struct_enum.in -> platt1, parent lattice A0
platt = readlatt(enumfile)
print 'Parent lattice';print platt
#enumerate lattice
os.chdir(maindir)
#os.system('enum.x %s' % enumfile)
##get scale from AFLOW
#scale = getscale(atomic,structchar)
#get scale from POSCAR of unit cell volume factor 1, in maindir
scale = getscalePOSCAR()
if scale < 0: #convert to simple length scale
    scale = (-scale/abs(det(platt)))**(1/3.0)
print 'Length scale for poscar', scale
#read struct list
structs = readstructs(structfile)

#multlist = [3,5]
#multlist = [2,4,6]
#multlist = [7,8]
#labeledstructs = [structchar + struct for struct in structs]
for struct in structs: #these are simply the numbers with no prefix
    for m in multlist:
        dir = structchar + struct + '_%i/' % m
        print;print dir
        if not os.path.isdir(finalDir+dir): os.system('mkdir %s' % finalDir+dir)              
        os.chdir(finalDir+dir)
        makestr2poscar(struct)
        if m == multlist[0] : #things to do only once per structure
            A = readlatt('POSCAR')
            [symops,nops] = getGroup(A)
            lattype = latticeType(nops)
            B = transpose(inv(A))
            [Bsymops,Bnops] = getGroup(A)
            L = getL(platt) #L is integer matrix for this lattice
            print 'Integer L matrix for lattice'; print L
            dL = abs(det(L))
            print 'det(L)', dL
            print 'detL*inv(L))'; print dL*inv(L)
            print 'Cubic mesh compatibility:'
            print 'detL*inv(L)*inv(platt)'; print dL*dot(inv(L),transpose(inv(platt)))
            print 'Fcc mesh compatibility for fcc lattice):'
            print 'detL*bccsq*trans(inv(L))'; print dL*dot(bccsq,transpose(inv(L)))
            
        scaleposcar(scale) #inserts length scale into appropriate spot.
        writefile([str(dL)],finalDir+dir+'detL')
        writefile([lattype],finalDir+dir+'lattype')
        os.system ('cp %s* %s' % (vaspinputdir,'.'))
        os.system ('rm slurm*')
        n = m*dL; shift = '0.00'       
        if meshtype == 'cub':
            M = n*transpose(inv(A))
            print 'M matrix for cubic mesh'; print M 
            print 'det M (kpoints in unreduced 1BZ)', det(M)            
            writekpts_cubic_n(n,shift)           
        if meshtype == 'fcc':
            M = trimSmall(n*dot(bccsq,transpose(inv(L))))
            print 'test K';print trimSmall(dot(B,inv(M)))
#            print 'test Bo';print transpose(inv(platt))
            print 'M matrix for fcc mesh'; print M
            print 'det M (kpoints in unreduced 1BZ)', det(M)
            
            writekpts_fcc_n(n,shift)
        writefile([str(abs(det(M)))],finalDir+dir+'detM')
        
#        writekpts_vasp_M(finalDir+dir,B,M,Bnops,Bsymops) #define K's explicitly!!!

        writejobfile(finalDir+dir)
        if n <= nmax:#(don't submit the biggest jobs
            ''
            
            subprocess.call(['sbatch', 'vaspjob']) #!!!!!!! Submit jobs
print 'done'
