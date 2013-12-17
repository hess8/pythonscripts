################# functions #######################
from numpy import array, cos, sin,arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import sign, matrix, transpose,rint,inner,multiply,size,argmin,round
from numpy import zeros 
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy import int as np_int
from numpy import float as np_float
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int
import time, os, subprocess, sys

#def structures_dict():
#    '''Dictionary for all 26 k-lattice cases in    
#            Setyawan, Wahyu; Curtarolo, Stefano (2010). "High-throughput electronic band structure calculations:   
#        {key:value}
#        {[number,'name',
#    '''
#    print
def isinteger(x):
    return np.equal(np.mod(x, 1), 0)

def isequal(x,y):
    eps = 1.0e-6
    return abs(x-y)<eps

def isreal(x):
    eps = 1.0e-6
    return abs(x.imag)<eps

def isindependent(vec1,vec2):  
    return not isequal(abs(cosvecs(vec1,vec2)),1)

def trimSmall(mat):
    low_values_indices = abs(mat) < 1.0e-6
    mat[low_values_indices] = 0.0
    return mat

def cosvecs(vec1,vec2):
    return dot(vec1,vec2)/norm(vec1)/norm(vec2)

class lattice(object): #reciprocal lattice
    def __init__(self):         
        self.vecs = []
        self.det = []
        self.Nmesh = []
        self.symops = []
        self.nops = []
        
def surfvol(vecs): 
    '''Surface/volume metric for a mesh, scaled to 1 for a cube'''
    vecs = transpose(vecs)
    u = norm(cross(vecs[0,:],vecs[1,:]))
    v = norm(cross(vecs[1,:],vecs[2,:]))
    w = norm(cross(vecs[2,:],vecs[0,:]))
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
                print 'sqrt%i ' % m, x
    return irrat

def abcalbega_latt(lvecs):
    '''finds a,b,c,alpha,beta,gamma from lattice vectors'''
    v0 = lvecs[0,:]
    v1 = lvecs[1,:]
    v2 = lvecs[2,:]
    a = norm(v0)
    b = norm(v1)
    c = norm(v2)
    alpha = arccos(dot(v1,v2)/b/c) * 180/pi
    beta = arccos(dot(v2,v0)/c/a) * 180/pi
    gamma = arccos(dot(v0,v1)/a/b) * 180/pi    
    return array([a,b,c,alpha,beta,gamma])

def readposcar(filename, path): 
    ''' Format is explicit lattice vectors, not a,b,c,alpha, beta, gamma'''
    file1 = open(path+filename,'r')
    poscar = file1.readlines()
    file1.close()
    descriptor = nstrip(poscar)[0]
    scale = float(nstrip(poscar)[1])
    if scale < 0:
        scale = (-scale)**(1/3.0)
    reallatt = zeros((3,3))
    reallatt[0,:] = array(poscar[2].split())
    reallatt[1,:] = array(poscar[3].split())
    reallatt[2,:] = array(poscar[4].split())
    reallatt = reallatt.astype(float)
    reallatt = scale*reallatt
    reciplatt = 2*pi*transpose(linalg.inv(reallatt))
    natoms = array(poscar[5].split(),dtype=int16)
    totatoms=sum(natoms)
    positions = zeros((totatoms,3),dtype=float)
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
    cryststruc = array(aflowin[istart+2].split(), dtype=float)
#    print scale
#    print cryststruc
    cryststruc[0:3] = scale*cryststruc[0:3]
    scale = 1.0 #since we put it into real lattice above
    print 'original a,b,c, alpha, beta, gamma'
    print cryststruc
    reallatt =  array(lattice_vecs(cryststruc)).astype(float)
    print 'reordered new lattice: a,b,c, alpha, beta, gamma'
    print abcalbega_latt(reallatt)

    print 'Lattice from aflow.in a,b,c alpha,beta,gamma > POSCAR0'
    print reallatt
    print
    reciplatt = 2*pi*transpose(linalg.inv(reallatt))
    print 'reciprocal lattice vectors'
    print reciplatt   
#    #test
#    print
#    print 'recip lattice traditional'
#    reciplatt2 = array(lattice_vecs(cryststruc)).astype(float)
#    a0 = reallatt[0,:]
#    a1 = reallatt[1,:]
#    a2 = reallatt[2,:]
#    vol = dot(a0,cross(a1,a2))
#    reciplatt2[0,:] = 2*pi*cross(a1,a2)/vol
#    reciplatt2[1,:] = 2*pi*cross(a2,a0)/vol
#    reciplatt2[2,:] = 2*pi*cross(a0,a1)/vol
#    print reciplatt2   
#    # end test
    natoms = array(aflowin[istart+3].split(),dtype=int16)
    totatoms=sum(natoms)
    positions = zeros((totatoms,3),dtype=float)
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
        [a,b,c,al,be,ga] = [a,c,b,al,ga,be] # e.g 3,2,1 -> 3,1,2
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
    lv[0,:] = [a,0,0]
    lv[1,:] = [b*cg,b*sg,0]  
    lv[2,0] = c*cb
    lv[2,1] = c/sg*(ca-cb*cg)
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
   
def svmesh(N,vecs):
    '''N: points desired.  vecs the lattice vectors as numpy array (reciprocal in our thinking)
    output:  n0, n1, n2, the number of divisions along each RLV for the mesh'''
    u = norm(cross(vecs[0,:],vecs[1,:]))
    v = norm(cross(vecs[1,:],vecs[2,:]))
    w = norm(cross(vecs[2,:],vecs[0,:]))
    n0 = N**(1/3.0) * w**(1/3.0) * u**(1/3.0) / v**(2/3.0)
    n1 = N**(1/3.0) * u**(1/3.0) * v**(1/3.0) / w**(2/3.0)
    n2 = N**(1/3.0) * v**(1/3.0) * w**(1/3.0) / u**(2/3.0)
    ns = [n0,n1,n2]
    print 'n0,n1,n2 from min s/v'
    print ns
    print

    p = n1/n0
    q = n2/n1
    r = n0/n2
    pqr = [p,q,r]
    print 'ratios n1/n0, n2/n1, n0/n2' 
    print pqr
    print
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
        r1i = argmax(abs(PQR)) #want index of the only nonzero element
        r1 = PQR[r1i]
        if r1 > 0: #CCW, so the higher index is greater m
            ms[r1i] = round(ns[r1i],0) #round smaller one
            ms[icy(r1i,1)] = ms[r1i] * r1 #form larger by relation
        else: #CW
            ms[icy(r1i,1)] = round(ns[icy(r1i,1)],0) #round smaller one
            ms[r1i] = ms[icy(r1i,1)] * abs(r1)      #form larger by relation
        ms[icy(r1i,-1)] = rint(ns[icy(r1i,-1)])
    if Nrels == 2 or Nrels == 3:
        #find the triangle side with no multiple relation.  Find the largest 
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
