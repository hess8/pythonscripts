################# functions #######################
import numpy as np   
import time, os, subprocess, sys

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
    print [a,b,c,al,be,ga]
    ca = np.cos(al/180*np.pi)
    cb = np.cos(be/180*np.pi)
    cg = np.cos(ga/180*np.pi)
    sa = np.sin(al/180*np.pi)
    sb = np.sin(be/180*np.pi)
    sg = np.sin(ga/180*np.pi) 
    lv = np.zeros((3,3))
    lv[0,:] = [a,0,0]
    lv[1,:] = [b*cg,b*sg,0]  
    lv[2,0] = c*cb
    lv[2,1] = c/sg*(ca-cb*cg)
    lv[2,2] = c/sg*np.sqrt(sg**2 - ca**2 - cb**2 + 2*ca*cb*cg)
    print 'sqrt', (sg**2 - ca**2 - cb**2 + 2*ca*cb*cg)
    return np.round(lv,14)       

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
   
def readposcar(path):    
    file1 = open(path+'/'+'POSCAR','r')
    poscar = file1.readlines()
    file1.close()
    natoms = np.sum(np.array(poscar[5].split(),dtype=np.int16))
    scale = float(poscar[1])
    if scale < 0:
        scale = np.abs(scale)**(1/3)
    reallatt = np.zeros((3,3))
    reallatt[0,:] = np.array(poscar[2].split())
    reallatt[1,:] = np.array(poscar[3].split())
    reallatt[2,:] = np.array(poscar[4].split())
    reallatt = scale*reallatt.astype(np.float)        
    reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
    return [natoms,reallatt,reciplatt]

def svmesh(N,vecs):
    '''N: points desired.  vecs the lattice vectors as numpy array (reciprocal in our thinking)
    output:  n0, n1, n2, the number of divisions along each RLV for the mesh'''
    u = np.linalg.norm(np.cross(vecs[0,:],vecs[1,:]))
    v = np.linalg.norm(np.cross(vecs[1,:],vecs[2,:]))
    w = np.linalg.norm(np.cross(vecs[2,:],vecs[0,:]))
    n0 = N**(1/3.0) * u**(1/3.0) * w**(1/3.0) / v**(2/3.0)
    n1 = N**(1/3.0) * u**(1/3.0) * v**(1/3.0) / w**(2/3.0)
    n2 = N**(1/3.0) * v**(1/3.0) * w**(1/3.0) / u**(2/3.0)
    ns = [n0,n1,n2]

    p = n1/n0
    q = n2/n1
    r = n0/n2
    pqr = [p,q,r]
    delta = 10**-6
    PQR = np.array([0,0,0])
    ms = np.array([0,0,0])
    '''   Define triangle as 
               m1
            P      R
        
        m2     Q     m3
    The integer relations (P,Q,R) 'point' to the larger integer.  If they are CCW, record +.  
    If CW, -    
           ''' 
    if abs(np.rint(p)-p)<delta:
        PQR[0] = np.rint(p)
    elif abs(np.rint(1/p)-(1/p))<delta:
        PQR[0] = -np.rint(1/p)
    if abs(np.rint(q)-q)<delta:
        PQR[1] = np.rint(q)
    elif abs(np.rint(1/q)-(1/q))<delta:
        PQR[1] = -np.rint(1/q)   
    if abs(np.rint(r)-r)<delta:
        PQR[2] = np.rint(r)
    elif abs(np.rint(1/r)-(1/r))<delta:
        PQR[2] = -np.rint(1/r)
    PQR = [int(PQR[j]) for j in [0,1,2]]
    Nrels = int(round(np.sum(np.abs(np.sign(PQR))),0)) #number of integer relations)
#    print 'Nrels', Nrels,PQR
#    print ns
#    print '0:1', ns[0]/ns[1], ns[1]/ns[0]
#    print '1:2', ns[1]/ns[2], ns[2]/ns[1]
#    print '2:0', ns[2]/ns[0], ns[0]/ns[2]
    #form final mesh m's
    if Nrels == 0:
        for i in [0,1,2]:
            ms[i] = np.rint(ns[i])
    if Nrels == 1:
        r1i = np.argmax(np.abs(PQR)) #want index of the only nonzero element
        r1 = PQR[r1i]
        if r1 > 0: #CCW, so the higher index is greater m
            ms[r1i] = round(ns[r1i],0) #round smaller one
            ms[icy(r1i,1)] = ms[r1i] * r1 #form larger by relation
        else: #CW
            ms[icy(r1i,1)] = round(ns[icy(r1i,1)],0) #round smaller one
            ms[r1i] = ms[icy(r1i,1)] * abs(r1)      #form larger by relation
        ms[icy(r1i,-1)] = np.rint(ns[icy(r1i,-1)])
    if Nrels == 2 or Nrels == 3:
        #find the triangle side with no integer relation.  Find the largest 
        nmax = np.max(np.abs(ns))
        imax = np.argmax(np.abs(ns)) #index of largest n.  Create m from this first
        r1 = abs(PQR[icy(imax,-1)])
        r2 = abs(PQR[imax])
        ms[imax] = r1 * r2 * int(np.rint(nmax/r1/r2))
        ms[icy(imax,-1)] = ms[imax]//r1
        ms[icy(imax,1)] = ms[imax]//r2
    return ms
                      
def getkpts_vasp(path):
    file1 = open(path+'/'+'KPOINTS','r')
    kpointsfile = file1.readlines()
    file1.close()
    if len(kpointsfile)>3:
        kpts_vasp = np.array(kpointsfile[3].split(),dtype=np.int32)
    else:
        kpts_vasp = []
    return kpts_vasp

def writekpts_vasp(dir, mesh):
    '''Write mesh m's to kpoints file, replacing previous mesh'''
    file1 = open(maindir+dir+'/'+kptsfile,'r')
    kpointsfile = file1.readlines()
    file1.close
    file2 = open(maindir+dir+'/'+kptsfile,'w')
    kpointsfile[0] = 'Mesh generated by s/v and integer relations method  \n'    
    kpointsfile[3] = '  '.join([str(mesh[i]) for i in [0,1,2]])+'\n'
    file2.writelines(kpointsfile) 
    file2.close()
    return 

def writejobfile(maindir,dir):
    '''read from a template in maindir, and put dir in job name'''
    file1 = open(maindir+'vaspjob','r')
    jobfile = file1.readlines()
    file1.close
    for i in range(len(jobfile)):
        jobfile[i]=jobfile[i].replace('myjob', dir)
    file2 = open(maindir+dir+'/'+'vaspjob','w')
    file2.writelines(jobfile) 
    file2.close()
    return 