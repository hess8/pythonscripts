'''Convention here is COLUMNS of matrices as vectors'''
################# functions #######################
from numpy import array, cos, sin,arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import sign, matrix, transpose,rint,inner,multiply,size,argmin,argmax,round,ceil
from numpy import zeros,nonzero,float64, sort, argsort,mod
fprec=float64
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy import int as np_int
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int
import time, os, subprocess, sys

utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') 
#had to copy and rename Gus's routine to the one below because ctypes could never find the one with the right name
getLatticePointGroup = utilslib.symmetry_module_mp_get_pointgroup_

def latticeType(nops):
    type = {2:'Triclinic', 4:'Monoclinic', 8:'Orthorhombic', 
            16:'Tetragonal', 12:'Trigonal', 24:'Hexagonal', 48:'Cubic'}
    return type[nops]

def unique_anorms(latt):
    '''Returns whether each vector in turn has a unique length among the three, '''
    a0 = norm(latt[:,0]); a1 = norm(latt[:,1]); a2 = norm(latt[:,2])
    return [not isequal(a0,a1) and  not isequal(a0,a2), not isequal(a1,a0) and not isequal(a1,a2), not isequal(a2,a0) and not isequal(a2,a1)] 

def searchsphere(aVecs):
   '''Decide how many lattice points to look in each direction to get all the
   points in a sphere that contains all of the longest _primitive_ vectors (from GLH)'''
   scale = 1  #1 is normal; others for testing
   eps  = 1e-6
   cell_volume = det(aVecs)
   max_norm = max(norm(aVecs[:,0]),norm(aVecs[:,1]),norm(aVecs[:,2]))
   n0 = scale*ceil(max_norm*norm(cross(aVecs[:,1],aVecs[:,2])/cell_volume)+eps)
   n1 = scale*ceil(max_norm*norm(cross(aVecs[:,2],aVecs[:,0])/cell_volume)+eps)
   n2 = scale*ceil(max_norm*norm(cross(aVecs[:,0],aVecs[:,1])/cell_volume)+eps)
   return [int(n0),int(n1),int(n2)]

def dNN(latt):
    '''Finds nearest neighbor distance'''
    [m0,m1,m2] = searchsphere(latt)
    dmin = 10000.0
    for i in range(-m0, m0):
        for j in range(-m1, m1):
                for k in range(-m2, m2):
                    d = i*latt[:,0] + j*latt[:,1] + k*latt[:,2]
                    dnorm = norm(d)
                    if dnorm < dmin and dnorm > 0.0:                    
                        dmin = dnorm
    return dmin

def lattvec_u(A,u):
    '''Find shortest lattice vector parallel to vector u.
    If v=|v|u is a lattice vector, made by multiplying unit vector by a real number |v|, then
    v = |v|u = m1a1 + m2a2 + m3a3 = A*col(m1 m2 m3), where col means column vector.  
    Then c === inv(A)u = col(m1/q m2/q m3/q).  We want to find the smallest set of m1, m2, m3 that satisfy this
    So divide c by the magnitude of the smallest of the nonzero elements of c (suppose it's the second).  
    Gives c2 = col(m1/m2 1 m3/m2).  We multiply this by the series m = 2, 3, 4, 5, until all elements are
    integers.  Then then c3  = m*c2 will be [m1 m2 m3] '''
    
#    print u
    c = trimSmall(dot(inv(A),u))
    c2 = c/min(abs(c[nonzero(c)]))
    c3 = c2
    mmax = 1000
    m = 1
#    print c3
    while not ( isinteger(c3[0]) and isinteger(c3[1]) and isinteger(c3[2]) ) and m < mmax :
#        print 'c3', c3
#        print 'mult',m,c3 
             
        m += 1
        c3 = m * c2

    if m < mmax: #success.
        v = dot(A,transpose(c3))
        return v
    else:
        return [0,0,0]

#    [m0,m1,m2] = searchsphere(latt)
#    scale = 2; [m0,m1,m2] = 2*array([m0,m1,m2])   #Needs 2x the value found in Gus' sphere routine
#    dmin = 10000.0
#    dvec = zeros((3),dtype = float)
#    lattvec = [0,0,0]
#    for i in range(-m0, m0):
#        for j in range(-m1, m1):
#                for k in range(-m2, m2):
#                    if i == j == k == 0: continue
#                    r = i*latt[:,0] + j*latt[:,1] + k*latt[:,2]
#                    d  = norm(r)
#                    ur = r/norm(r)
#                    if 0.00 < d < dmin and isequal(dot(transpose(ur),u),1):                    
#                        dmin = d
#                        lattvec = r             
    return lattvec

#def three_perp(latt,lattype):
#    '''Find three of the shortest lattice vectors that are perpendicular to each other'''
#    [m0,m1,m2] = searchsphere(latt)
##    print [m0,m1,m2]
#    vecs = zeros(((2*m0+1) * (2*m1+1) * (2*m2+1),3), dtype=float)
#    norms = zeros((2*m0+1) * (2*m1+1) * (2*m2+1), dtype=float) 
#    index = 0
#    #find all vectors in this range
#    for i in range(-m0, m0+1):
#        for j in range(-m1, m1+1):
#                for k in range(-m2, m2+1):
##                    print i,k,j, index
#                    vecs[index,:] = i*latt[:,0] + j*latt[:,1] + k*latt[:,2]
##                    print  vecs[index,:],norm(vecs[index,:])
#                    norms[index]  = norm(vecs[index,:])
#                    index += 1
#    #Find indices for vectors sorted by increasing norm
#    vecs = trimSmall(vecs)
#    sortvecs = argsort(norms)
##    for i in sortvecs:
##        print vecs[i,:]
##    print sortvecs
##    print sort(norms)
##    for i in sortvecs:
##        print vecs[i,:],
##        print norms[i]
#    #Find the first set of three shortest that are perpendicular and meet symmetry requirements
#    for i in sortvecs[1:]: #Don't include first which i zero
#        for j in sortvecs[1:]:
#            if arenormal(vecs[i,:],vecs[j,:]):
##                print i,j       
#                for k in sortvecs[1:]:
##                    if i != j and i!=k and k != j:
##                        print;print k, cosvecs(vecs[k,:],vecs[j,:]) , cosvecs(vecs[k,:],vecs[i,:])
#                    if arenormal(vecs[k,:],vecs[j,:]) and arenormal(vecs[i,:],vecs[k,:]):
#                        S = transpose(array([vecs[i,:],vecs[j,:],vecs[k,:]]))
#                        print 'found 3 vectors of any type',i,j,k
#                        if lattype == 'Cubic' and unique_anorms(S).count(True) == 0:
#                            return trimSmall(S) 
#                        if lattype == 'Tetrahedral' and unique_anorms(S).count(True) == 1:
#                            return trimSmall(S) 
#                        if lattype == 'Orthorhombic' and unique_anorms(S).count(True) == 3:
#                            return trimSmall(S)
##                    print vecs[i,:];print vecs[j,:]; print vecs[k,:]
##                    sys.exit('stop in 3perp')   
#    sys.exit('Could not find 3 perp vectors')

def packingFraction(latt):
#    print 'dNN', dNN(latt)
    vol = abs(det(latt))
    vsphere = 4/3.0*pi*(dNN(latt)/2)**3.0  
    return round(vsphere/vol,4)

def isinteger(x):
    return isequal(abs(rint(x)-x), 0)

def isequal(x,y):
    eps = 5.0e-5
    return abs(x-y)<eps

def isreal(x):
    eps = 1.0e-6
    return abs(x.imag)<eps

def arenormal(v1,v2):
    eps = 5.0e-5
    norm1 = norm(v1); norm2 = norm(v2);
    if norm1*norm2>eps:
#        print cosvecs(v1,v2)
        return abs(dot(v1,v2))<eps*norm1*norm2
    else:
        return False

def isindependent(vec1,vec2):  
    return not isequal(abs(cosvecs(vec1,vec2)),1)

def trimSmall(list_mat):
    low_values_indices = abs(list_mat) < 1.0e-5
    list_mat[low_values_indices] = 0.0
    return list_mat

def cosvecs(vec1,vec2):
    return dot(vec1,vec2)/norm(vec1)/norm(vec2)

class lattice(object): #reciprocal lattice
    def __init__(self):         
        self.vecs = []
        self.det = []
        self.Nmesh = []
        self.symops = []
        self.msymops = []
        self.nops = []
        self.lattype = []
        
def surfvol(vecs): 
    '''Surface/volume metric for a mesh, scaled to 1 for a cube'''
    vecs = transpose(vecs)
    u = norm(cross(vecs[:,0],vecs[:,1]))
    v = norm(cross(vecs[:,1],vecs[:,2]))
    w = norm(cross(vecs[:,2],vecs[:,0]))
    surfvol = 2*(u + v + w)/6/(abs(det(vecs)))**(2/3.0)
    return abs(surfvol)

def orthdef(latt): 
    '''Orthogonality defect''' #columns as vectors
    od = norm(latt[:,0])*norm(latt[:,1])*norm(latt[:,2])/det(latt)
    return abs(od)

def nstrip(list):
#    '''Strips off /n'''
    import string
    list2 = []
    for string1 in list:   
        string2 = string1.strip("\n")
        list2.append(string2)
    return list2

def intcheck(x):
    delta = 10**-6
    if abs(rint(x)-x)<delta:
        return True
    else:
        return False

def irratcheck(ratios,mlist):
    '''Checks to see if ratios are a multiple of a root, or a multiple of 1/root'''
    irrat = ''
    for m in mlist:                 
        sqr = sqrt(m)
        for i,x in enumerate(ratios):
            if intcheck(x*sqr) or intcheck(x/sqr) or intcheck(sqr/x) or intcheck((1/x*sqr)):
                irrat = irrat+'sqrt%i ' % m
#                print 'sqrt%i ' % m, x
    return irrat

def abcalbega_latt(lvecs):
    '''finds a,b,c,alpha,beta,gamma from lattice vectors'''
    v0 = lvecs[:,0]
    v1 = lvecs[:,1]
    v2 = lvecs[:,2]
    a = norm(v0)
    b = norm(v1)
    c = norm(v2)
    alpha = arccos(dot(v1,v2)/b/c) * 180/pi
    beta = arccos(dot(v2,v0)/c/a) * 180/pi
    gamma = arccos(dot(v0,v1)/a/b) * 180/pi    
    return array([a,b,c,alpha,beta,gamma])

def readposcar(filename, path): 
    ''' Format is explicit lattice vectors, not a,b,c,alpha, beta, gamma. 
    Saves vectors as columns, not rows which POSCAR uses'''
    file1 = open(path+filename,'r')
    poscar = file1.readlines()
    file1.close()
    descriptor = nstrip(poscar)[0]
    scale = float(nstrip(poscar)[1])
    reallatt = zeros((3,3),dtype=fprec)
    reallatt[:,0] = array(poscar[2].split())
    reallatt[:,1] = array(poscar[3].split())
    reallatt[:,2] = array(poscar[4].split())
#    reallatt = reallatt.astype(fprec)
#    print reallatt
    if scale < 0:
        vol = det(reallatt)
        scale = (-scale/vol)**(1/3.0)  
    reallatt = scale*reallatt
#    print reallatt    
    scale = 1.0
    reciplatt = 2*pi*transpose(inv(reallatt))
    natoms = array(poscar[5].split(),dtype=int)
    totatoms=sum(natoms)
    positions = zeros((totatoms,3),dtype=fprec)
    postype = poscar[6] #Direct or Cartesian
    whichatom = 0
    for natom in natoms:
        for i in range(natom):
            for k in [0,1,2]:
                positions[whichatom,k] = float(poscar[7+whichatom].split()[k])
            whichatom += 1
    totatoms=sum(natoms)
    return [descriptor, scale, reallatt, reciplatt, natoms, postype, positions]


def aflow2poscar(path): 
    '''DON"T USE THIS>>>Need to use (-scale/volume)**(1/3.0)'''
    file1 = open(path+'aflow.in','r')
    aflowin = file1.readlines()
    file1.close()
    for i,line in enumerate(aflowin):
        if 'VASP_POSCAR_MODE_EXPLICIT' in line:
            istart = i+1
            break #take only first instance (should be only one)
    descriptor = nstrip(aflowin)[istart]
    scale = float(nstrip(aflowin)[istart+1])
    if scale < 0:
        scale = (-scale)**(1/3.0)
    cryststruc = array(aflowin[istart+2].split(), dtype=fprec)
#    print scale
#    print cryststruc
    cryststruc[0:3] = scale*cryststruc[0:3]
    scale = 1.0 #since we put it into real lattice above
    print 'original a,b,c, alpha, beta, gamma'
    print cryststruc
    reallatt =  array(lattice_vecs(cryststruc)).astype(fprec)
    print 'reordered new lattice: a,b,c, alpha, beta, gamma'
    print abcalbega_latt(reallatt)

    print 'Lattice from aflow.in a,b,c alpha,beta,gamma > POSCAR0'
    print reallatt
    print
    reciplatt = 2*pi*transpose(inv(reallatt))
    print 'reciprocal lattice vectors'
    print reciplatt   
#    #test
#    prints
#    print 'recip lattice traditional'
#    reciplatt2 = array(lattice_vecs(cryststruc)).astype(fprec)
#    a0 = reallatt[:,0]
#    a1 = reallatt[:,1]
#    a2 = reallatt[:,2]
#    vol = dot(a0,cross(a1,a2))
#    reciplatt2[:,0] = 2*pi*cross(a1,a2)/vol
#    reciplatt2[:,1] = 2*pi*cross(a2,a0)/vol
#    reciplatt2[:,2] = 2*pi*cross(a0,a1)/vol
#    print reciplatt2   
#    # end test
    natoms = array(aflowin[istart+3].split(),dtype=int)
    totatoms=sum(natoms)
    positions = zeros((totatoms,3),dtype=fprec)
    postype = aflowin[istart+4] #Direct or Cartesian
    where = 0
    for natom in natoms:
        for i in range(natom):
            for k in [0,1,2]:
                positions[where,k] = float(aflowin[istart+5+where].split()[k])
            where += 1
    create_poscar('POSCAR',descriptor+'From aflow.in BCH',scale,reallatt,natoms,postype,positions,path)
    totatoms=sum(natoms)
    return totatoms

def create_poscar(filename,descriptor, scale, latticevecs, natoms, type_pos, positions, path):
    poscar = open(path+filename,'w')
    poscar.write(descriptor+'\n')
    poscar.write(str(scale)+'\n')
    for i in [0,1,2]:
        poscar.write('%20.15f %20.15f %20.15f \n' % (latticevecs[i,0], latticevecs[i,1], latticevecs[i,2]))         
    for i in natoms:
        poscar.write(str(i)+'    ')
    poscar.write('\n')
    poscar.write(type_pos)
    where = 0
    for natom in natoms:
        for i in range(natom):
            poscar.write('%20.15f %20.15f %20.15f \n' % (positions[where,0],positions[where,1],positions[where,2]))
            where += 1
    poscar.close()

def lattice_vecs(cryststruc):
    ''' Make lattice vectors from triclinic method of 
        Setyawan, Wahyu; Curtarolo, Stefano (2010). "High-throughput electronic band structure calculations:   
    see A.14.
    al, be, ga, are alpha, beta, gamma angles
    Stefano's method reorders so that a<b<c .  Then the angles are changed:
    a-b:gamma  b-c: alph  c-a: beta
    '''
    [a,b,c,al,be,ga] = cryststruc
    if b>c: #switch c,b, gamma and beta 
        [a,b,c,al,be,ga] = [a,c,b,al,ga,be] # e.g 3,1,2 -> 3,1,2
    if a>b: #switch a,b, alph and beta
        [a,b,c,al,be,ga] = [b,a,c,be,al,ga] # e.g 3,1,2 -> 1,3,2
    if b>c: #switch c,b, gamma and beta
        [a,b,c,al,be,ga] = [a,c,b,al,ga,be] # e.g 1,3,2 -> 1,2,3
#    print [a,b,c,al,be,ga]
    ca = cos(al/180.0*pi)
    cb = cos(be/180.0*pi)
    cg = cos(ga/180.0*pi)
    sa = sin(al/180.0*pi)
    sb = sin(be/180.0*pi)
    sg = sin(ga/180.0*pi) 
    lv = zeros((3,3))
    lv[:,0] = [a,0,0]
    lv[:,1] = [b*cg,b*sg,0]  
    lv[0,2] = c*cb
    lv[1,2] = c/sg*(ca-cb*cg)
    lv[2,2] = c/sg*sqrt(sg**2 - ca**2 - cb**2 + 2*ca*cb*cg)
    return round(lv,14)       

def icy(i,change): #for cycling indices 0,1,2
    i = i+change
    if i>2:
        i=0
    if i<0:
        i=2
    return i

def regpy_nocase(str,path):
    import re
    file1 = open(path,'r')
    lines = file1.readlines()
    file1.close()
    for line in lines:
        if re.search( str, line,  re.M|re.I):
            print line

def svmesh1freedir(N,vecs):
    '''N: points desired, when only one direction is free of contraints on its magnitude vs other directions.
    Put this direction first in vecs (vecs[:,0])!!  Vecs the lattice vectors as numpy array (reciprocal in our thinking)
    output:  n0, n1=n2, the number of divisions along each RLV for the mesh. '''
    u = norm(cross(vecs[:,0],vecs[:,1])) 
    v = norm(cross(vecs[:,1],vecs[:,2]))
    w = norm(cross(vecs[:,2],vecs[:,0]))
    n0 = rint(N**(1/3.0) * (2*(u + w)/v)**(2/3.0))
    n1 = rint(sqrt(N/n0))
    return [n0,n1] 
   
def svmesh(N,vecs):
    '''N: points desired, when all three directions are free of constraints on its magnitude vs other directions.  
    Vecs the lattice vectors as numpy array (reciprocal here)
    output:  n0, n1, n2, the number of divisions along each RLV for the mesh'''
#    print 'vecs in svmesh', vecs
    u = norm(cross(vecs[:,0],vecs[:,1]))
    v = norm(cross(vecs[:,1],vecs[:,2]))
    w = norm(cross(vecs[:,2],vecs[:,0]))
    n0 = N**(1/3.0) * w**(1/3.0) * u**(1/3.0) / v**(2/3.0)
    n1 = N**(1/3.0) * u**(1/3.0) * v**(1/3.0) / w**(2/3.0)
    n2 = N**(1/3.0) * v**(1/3.0) * w**(1/3.0) / u**(2/3.0)
    ns = array([n0,n1,n2])
#    print 'n0,n1,n2 from min s/v'
#    print ns
#    print

    p = n1/n0
    q = n2/n1
    r = n0/n2
    pqr = [p,q,r]
#    print 'ratios n1/n0, n2/n1, n0/n2' 
#    print pqr
#    print
    primes = [2,3,5,7,11,13,17,19,23,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]
    irrat = irratcheck(pqr,primes)

    PQR = array([0,0,0])
    ms = array([0,0,0])
    '''   Define triangle as 
               m1
            P      R
        
        m2     Q     m3
    The multiple relations (P,Q,R) 'point' to the larger number.  If they are CCW, record +.  
    If CW, -    
           ''' 
    if intcheck(p):
        PQR[0] = rint(p)
    elif intcheck(1/p):
        PQR[0] = -rint(1/p)
    if intcheck(q):
        PQR[1] = rint(q)
    elif intcheck(1/q):
        PQR[1] = -rint(1/q)        
    if intcheck(r):
        PQR[2] = rint(r)
    elif intcheck(1/r):
        PQR[2] = -rint(1/r)        
    PQR = [int(PQR[j]) for j in [0,1,2]]
    Nrels = int(round(sum(abs(sign(PQR))),0)) #number of multiple relations)
#    print 'Nrels', Nrels,PQR
#    print ns
#    print '0:1', ns[0]/ns[1], ns[1]/ns[0]
#    print '1:2', ns[1]/ns[2], ns[2]/ns[1]
#    print '2:0', ns[2]/ns[0], ns[0]/ns[2]
    #form final mesh m's
    if Nrels == 0:
        for i in [0,1,2]:
            ms[i] = rint(ns[i])
    if Nrels == 1:
        r1i = argmax(abs(array(PQR))) #want index of the only nonzero element
        r1 = PQR[r1i]
        if r1 > 0: #CCW, so the higher index is greater m
            ms[r1i] = round(ns[r1i],0) #round smaller one
            ms[icy(r1i,1)] = ms[r1i] * r1 #form larger by relation
        else: #CW
            ms[icy(r1i,1)] = round(ns[icy(r1i,1)],0) #round smaller one
            ms[r1i] = ms[icy(r1i,1)] * abs(r1)      #form larger by relation
        ms[icy(r1i,-1)] = rint(ns[icy(r1i,-1)])
    if Nrels == 2:#relations either point away from the smallest or toward the largest n
        zindex = PQR.index(0)#find the triangle side with no multiple relation.
        mpivot = icy(zindex,-1)#the vertex opposite a side is -1 away
        if mpivot == argmax(abs(ns)):
            nmax = max(abs(ns))
            imax = mpivot #index of largest n.  Create m from this first
            r1 = abs(PQR[icy(imax,-1)])
            r2 = abs(PQR[imax])
            ms[imax] = r1 * r2 * int(rint(nmax/r1/r2))
            ms[icy(imax,-1)] = ms[imax]/r1
            ms[icy(imax,1)] = ms[imax]/r2  
        if mpivot == argmin(abs(ns)):
            nmin = min(abs(ns))
            imin = mpivot #index of largest n.  Create m from this first
            r1 = abs(PQR[icy(imin,-1)])
            r2 = abs(PQR[imin])
            p
            ms[imin] = int(rint(nmin))
            ms[icy(imin,-1)] = ms[imin]*r1
            ms[icy(imin,1)] = ms[imin]*r2                  
    if Nrels == 3: #either method used in 2 works
        nmax = max(abs(ns))
        imax = argmax(abs(ns)) #index of largest n.  Create m from this first
        r1 = abs(PQR[icy(imax,-1)])
        r2 = abs(PQR[imax])
        ms[imax] = r1 * r2 * int(rint(nmax/r1/r2))
        ms[icy(imax,-1)] = ms[imax]/r1
        ms[icy(imax,1)] = ms[imax]/r2
    return [ms,irrat]
                      
def getkpts_vasp(path):
    file1 = open(path+'/'+'KPOINTS','r')
    kpointsfile = file1.readlines()
    file1.close()
    if len(kpointsfile)>3:
        kpts_vasp = array(kpointsfile[3].split(),dtype=int32)
    else:
        kpts_vasp = []
    return kpts_vasp

def writekpts_vasp(maindir,dir,kptsfile,Nkppra):
    '''Write mesh m's to kpoints file, replacing previous mesh'''
    file1 = open(maindir+dir+'/'+kptsfile,'r')
    kpointsfile = file1.readlines()
    file1.close
    file2 = open(maindir+dir+'/'+kptsfile,'w')
#    kpointsfile[0] = kpointsfile[0].strip('\n') + ' BCH\n' 
    kpointsfile[0] = 'Generated by BCH using Nkppra = %s\n' %Nkppra
    kpointsfile[2] = 'Minkowski Monkhorst-Pack\n' 
    kpointsfile[3] = str(Nkppra)+'\n'
    file2.writelines(kpointsfile) 
    file2.close()
    return 

def writejobfile(maindir,struct,jobfile,nameadd,vaspexec):
    '''read from a template one level up from maindir, and put dir in job name ('myjob' is replaced).  
    Choose executable label (should be linked in bin file)'''
    curdir = os.getcwd()
    os.chdir(maindir)
    file1 = open('../'+jobfile,'r')
    jobfile = file1.readlines()
    file1.close
    for i in range(len(jobfile)):
        jobfile[i]=jobfile[i].replace('myjob', struct+nameadd)
        if '>' in jobfile[i]:
            jobfile[i]= vaspexec + ' > vasp.out'
    file2 = open(maindir+dir+'vaspjob','w')
    file2.writelines(jobfile) 
    file2.close()
    os.chdir(curdir)
    return 

def checkq(user):
    ####### This "if" block is for python versions lower than 2.7. Needed for subprocess.check_output. 
    if "check_output" not in dir( subprocess ): # duck punch it in!
        def f(*popenargs, **kwargs):
            if 'stdout' in kwargs:
                raise ValueError('stdout argument not allowed, it will be overridden.')
            process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
            output, unused_err = process.communicate()
            retcode = process.poll()
            if retcode:
                cmd = kwargs.get("args")
                if cmd is None:
                    cmd = popenargs[0]
                raise CalledProcessError(retcode, cmd)
            return output
        subprocess.check_output = f
    ####### above needed for python versions lower than 2.7.
################# script #######################
    maxDays = 4 #days to run this script
    waitMin = 0.2 #minutes between checking         
    #os.system('rm slurm-*.out') 
    starttime = time.time()    
    while time.time()-starttime < maxDays*3600*24: 
        pending = subprocess.check_output(['squeue','-u',user,'--state=PENDING'])
        running = subprocess.check_output(['squeue','-u',user,'--state=RUNNING'])    
    #    pending = subprocess.check_output(['squeue','-u',user,'-n',jobname,'--state=PENDING'])
    #    running = subprocess.check_output(['squeue','-u',user,'-n',jobname,'--state=RUNNING'])
    #    locked = subprocess.check_output(['find','-name','LOCK','|','wc','-l']) #gives error
        print pending
        pending = pending.splitlines()
        running = running.splitlines() 
        print 'Jobs that are pending: ',len(pending)-1 #one line is header     
        print 'Jobs that are running: ',len(running)-1 #one line is header 
        print 'Will check again in %s min' % waitMin
    #    print 'Locked files:' , locked   
        print
        time.sleep(waitMin*60) #seconds between checks
        
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
    a = zeros((3,3,nops),dtype=fprec)
    ielement = 0
    for i in range(3):
        for j in range(3):
            for k in range(nops):                          
                a[i,j,k] = OUT[ielement]
                ielement += 1                 
    return a

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
    symopsB = trimSmall(unload_ctypes_3x3xN_double(opsOUT,nops))
    return [symopsB,nops]

def intsymops(A):
    '''finds integer symmetry operations in the basis of the vectors of lattice A, 
    given the cartesian operators A.symop'''
    mops = zeros((3,3,A.nops),dtype = float)
    for i in range(A.nops):
        mops[:,:,i] = dot(dot(inv(A.vecs),A.symops[:,:,i]),A.vecs)
    return trimSmall(mops)  
      
def checksymmetry(latt,parentlatt):
    '''check that the lattice obeys all symmetry operations of a parent lattice:  R.latt.inv(R) will give an integer matrix'''
    for iop in range(parentlatt.nops):
        lmat = array(latt)
        if det(lmat) == 0:
            print 'Determinant zero'
            print lmat    
            return False
        mmat = trimSmall(dot(dot(inv(lmat),parentlatt.symops[:,:,iop]),lmat))
#        print 'mmat', iop
#        print trimSmall(mmat)
        for i in range(3):
            for j in range(3):
                if abs(rint(mmat[i,j])-mmat[i,j])>1.0e-4:
#                    print mmat
#                    print iop, 'Symmetry failed for mmat[i,j]',mmat[i,j]
#                    print 'Cartesian operator' 
#                    print parentlatt.symops[:,:,iop] 
#                    print 'Cartesian Lattice'
#                    print lmat
#                    print 'Noninteger operator in m space' 
#                    print mmat                                        
                    return False #jumps out of subroutine
    return True #passes test

def symmetryError(latt,parentlatt):
    '''check that the lattice obeys all symmetry operations of a parent lattice:  R.latt.inv(R) will give an integer matrix'''
    symmerr = 0.0
    for iop in range(parentlatt.nops):
        lmat = array(latt)
#        if det(lmat) == 0:
#            print 'Determinant zero'
#            print lmat    
#            return symmerror
        mmat = trimSmall(dot(dot(inv(lmat),parentlatt.symops[:,:,iop]),lmat))
#        print 'mmat', iop
#        print trimSmall(mmat)
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
        if operr > 2.0e-4:
            symmerr += operr
#            print 'Noninteger operator in superlattice for operation %s, with error %s.' % (iop,str(operr))
#    if operr < 2.0e-4:
#        return 0.0
#    else:
    return symmerr

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
   
def MT2mesh(MT,B,A):
    '''tests which directions are free to adjust, and chooses best mesh'''
    Nmesh = B.Nmesh
    testi = 7
    freedir = []
    M = transpose(MT)
#    print 'M';print M
#    print 'inv(M)';print inv(M)
    Q = dot(B.vecs,inv(M))    
    print 'starting mesh Q'; print trimSmall(Q)
    if checksymmetry(Q,B): 
        print 'Pass symmetry check' 
    else: 
        print'Failed symmetry check'
    for i in range(3):
        Qtest = dot(B.vecs,inv(M))
        Qtest[:,i] = Qtest[:,i]/testi
        if checksymmetry(Qtest,B): freedir.append(i)
    print 'free directions', freedir
    if len(freedir) == 0:
        N2 = rint((Nmesh/abs(det(M)))**(1/3.0))
        ms = [N2,N2,N2]    
    elif len(freedir) == 1:
        #order so free dir is first in vecs
        freeindex = freedir[0]
        otherindices = nonzero(array([0,1,2])-freeindex)
        vecs = zeros((3,3),dtype = fprec)
        vecs[:,0] = Q[:,freeindex]
        vecs[:,1] = Q[:,otherindices[0][0]]
        vecs[:,2] = Q[:,otherindices[0][1]]
        [n0,n1]= svmesh1freedir(Nmesh/abs(det(M)),vecs)
        print [n0,n1]
        ms = [0,0,0]
        ms[freeindex] = n0
        ms[otherindices[0][0]] = n1
        ms[otherindices[0][1]] = n1
    elif len(freedir) == 3:   
        mesh = svmesh(Nmesh/abs(det(M)),Q)
        ms = mesh[0]
    else: #should not occur
        test2free(M,B,A)
        sys.exit('Error in MT2mesh; freedir has 2 elements, but not 3')
    print 'mesh integers', ms      
    Q[:,0] = Q[:,0]/ms[0]
    Q[:,1] = Q[:,1]/ms[1]
    Q[:,2] = Q[:,2]/ms[2]
    return Q
#    if checksymmetry(Q,B):

def test2free(M,B,A):
    '''just for testing this strange case'''
    [m1,m2,m3] = [1,1,4]
    testi =7
    freedir = []
    
    Qtest = dot(B.vecs,inv(M))
    Qtest[:,0] = Qtest[:,0]/m1
    Qtest[:,1] = Qtest[:,1]/m2
    Qtest[:,2] = Qtest[:,2]/m3
    print 'Q symmetry test for B:', checksymmetry(Qtest,B) 
 
    Stest1 = dot(A.vecs,transpose(M))
    print; print 'starting superlattice S'; print trimSmall(Stest1)
    print 'Initial S symmetry test for A:', checksymmetry(Stest1,A)
    for i in range(3):
        Stest2 = dot(A.vecs,transpose(M))
        Stest2[:,i] = Stest2[:,i]*testi
        if checksymmetry(Stest2,A): freedir.append(i)
    print 'free directions for A', freedir
    Stest = dot(A.vecs,transpose(M))    
    Stest[:,0] = Stest[:,0]*m1
    Stest[:,1] = Stest[:,1]*m2
    Stest[:,2] = Stest[:,2]*m3
    print 'Final S'; print trimSmall(Stest)
    print 'Final S symmetry test for A:', checksymmetry(Stest,A)  
#        if not checksymmetry(Qtest,B):
#        print 'B lattice transpose';print 10*transpose(B.vecs)
#        print 'Q mesh transpose';print 10*transpose(Qtest)
#    if not checksymmetry(Stest,A):
#    print 'A lattice transpose for VMD';print trimSmall(transpose(A.vecs))
#    print 'S transpose for VMD';print trimSmall(transpose(Stest))
#    
#    print 'S transpose lattice after x rotation';print trimSmall(transpose(dot(A.symops[:,:,2],Stest)))
    mmat = trimSmall(dot(dot(inv(Stest),A.symops[:,:,2]),Stest))
    print 'mmat'; print mmat
    print 'RS,Sm';print trimSmall(dot(A.symops[:,:,2],Stest));print trimSmall(dot(Stest,mmat))
    mmat1 = trimSmall(dot(dot(inv(Stest1),A.symops[:,:,2]),Stest1))
    print 'mmat1'; print mmat1
    print 'RSinitial,Sinitial*m'; print trimSmall(dot(A.symops[:,:,2],Stest1));print trimSmall(dot(Stest1,mmat1))

#    print 'm from inverse':  