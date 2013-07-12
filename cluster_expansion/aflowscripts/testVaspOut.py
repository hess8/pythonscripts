#!/usr/bin/python
import time, os, subprocess
'''For each dir in jobs2run: replaces kpoints file with correct mesh for POSCAR,
reads a jobfiles from the maindir,writes the structure number to the job name, and submits a vasp job
'''
#!/usr/bin/python
''' tests. '''
    
import sys,os
import numpy as np
################# functions #######################
def icy(i,change): #for cycling indices 0,1,2
    i = i+change
    if i>2:
        i=0
    if i<0:
        i=2
    return i

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
    print 'Nrels', Nrels,PQR
    print ns
    
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
    print ms, 'ms'
    return ms
                      
def getkpts_vasp(dir):
    file1 = open(maindir+dir+'/'+kptsfile,'r')
    kpointsfile = file1.readlines()
    file1.close()
    kpts_vasp = np.array(kpointsfile[3].split(),dtype=np.int16)
    return kpts_vasp

def writekpts_vasp(dir, mesh):
    '''Write mesh m's to kpoints file, replacing previous mesh'''
    file1 = open(maindir+dir+'/'+kptsfile,'r')
    kpointsfile = file1.readlines()
    file1.close
    file2 = open(maindir+dir+'/'+kptsfile,'w')
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


################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAktest2corr/AlIr/'
testfile = 'POSCAR'
kptsfile = 'KPOINTS'
Nkppra = 10000

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
        print
        print dir
        file1 = open(maindir+dir+'/'+testfile,'r')
        poscar = file1.readlines()
        file1.close()
        if len(poscar) >0:
            scale = np.sum(np.array(float(poscar[1])))
            N = np.rint(Nkppra/np.sum(np.array(poscar[5].split(),dtype=np.int16))).astype(int) # number of kpts desired
            reallatt[0,:] = np.array(poscar[2].split())
            reallatt[1,:] = np.array(poscar[3].split())
            reallatt[2,:] = np.array(poscar[4].split())
            reallatt = scale*reallatt.astype(np.float)        
            reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
            mesh_ns = svmesh(N,reciplatt)
            writekpts_vasp(dir,mesh_ns) #correct kmesh
            writejobfile(maindir,dir)
            os.chdir(dir)
            subprocess.call(['rm', 'slurm*'])
            subprocess.call(['sbatch', 'vaspjob'])
            os.chdir(maindir)
            
            # submit vasp job                      
print 'Done'



