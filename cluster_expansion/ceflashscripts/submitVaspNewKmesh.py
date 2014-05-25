#!/usr/bin/python
import time, os, subprocess
'''For each dir in jobs2run: copies POSCAR.orig to POSCAR, replaces kpoints file with correct mesh for POSCAR,
reads a jobfiles from the maindir,writes the structure number to the job name, and submits a vasp job
'''
#!/usr/bin/python
''' tests. '''
    
import sys,os
import numpy as np
#from kmeshroutines.py import *
from kmeshroutines.py import getkpts_vasp
################# functions #######################
#def icy(i,change): #for cycling indices 0,1,2
#    i = i+change
#    if i>2:
#        i=0
#    if i<0:
#        i=2
#    return i
#
#def svmesh(N,vecs):
#    '''N: points desired.  vecs the lattice vectors as numpy array (reciprocal in our thinking)
#    output:  n0, n1, n2, the number of divisions along each RLV for the mesh'''
#    u = np.linalg.norm(np.cross(vecs[0,:],vecs[1,:]))
#    v = np.linalg.norm(np.cross(vecs[1,:],vecs[2,:]))
#    w = np.linalg.norm(np.cross(vecs[2,:],vecs[0,:]))
#    n0 = N**(1/3.0) * u**(1/3.0) * w**(1/3.0) / v**(2/3.0)
#    n1 = N**(1/3.0) * u**(1/3.0) * v**(1/3.0) / w**(2/3.0)
#    n2 = N**(1/3.0) * v**(1/3.0) * w**(1/3.0) / u**(2/3.0)
#    ns = [n0,n1,n2]
#
#    p = n1/n0
#    q = n2/n1
#    r = n0/n2
#    pqr = [p,q,r]
#    delta = 10**-6
#    PQR = np.array([0,0,0])
#    ms = np.array([0,0,0])
#    '''   Define triangle as 
#               m1
#            P      R
#        
#        m2     Q     m3
#    The integer relations (P,Q,R) 'point' to the larger integer.  If they are CCW, record +.  
#    If CW, -    
#           ''' 
#    if abs(np.rint(p)-p)<delta:
#        PQR[0] = np.rint(p)
#    elif abs(np.rint(1/p)-(1/p))<delta:
#        PQR[0] = -np.rint(1/p)
#    if abs(np.rint(q)-q)<delta:
#        PQR[1] = np.rint(q)
#    elif abs(np.rint(1/q)-(1/q))<delta:
#        PQR[1] = -np.rint(1/q)   
#    if abs(np.rint(r)-r)<delta:
#        PQR[2] = np.rint(r)
#    elif abs(np.rint(1/r)-(1/r))<delta:
#        PQR[2] = -np.rint(1/r)
#    PQR = [int(PQR[j]) for j in [0,1,2]]
#    Nrels = int(round(np.sum(np.abs(np.sign(PQR))),0)) #number of integer relations)
#    print 'Nrels', Nrels,PQR
#    print ns
#    
#    #form final mesh m's
#    if Nrels == 0:
#        for i in [0,1,2]:
#            ms[i] = np.rint(ns[i])
#    if Nrels == 1:
#        r1i = np.argmax(np.abs(PQR)) #want index of the only nonzero element
#        r1 = PQR[r1i]
#        if r1 > 0: #CCW, so the higher index is greater m
#            ms[r1i] = round(ns[r1i],0) #round smaller one
#            ms[icy(r1i,1)] = ms[r1i] * r1 #form larger by relation
#        else: #CW
#            ms[icy(r1i,1)] = round(ns[icy(r1i,1)],0) #round smaller one
#            ms[r1i] = ms[icy(r1i,1)] * abs(r1)      #form larger by relation
#        ms[icy(r1i,-1)] = np.rint(ns[icy(r1i,-1)])
#    if Nrels == 2 or Nrels == 3:
#        #find the triangle side with no integer relation.  Find the largest 
#        nmax = np.max(np.abs(ns))
#        imax = np.argmax(np.abs(ns)) #index of largest n.  Create m from this first
#        r1 = abs(PQR[icy(imax,-1)])
#        r2 = abs(PQR[imax])
#        ms[imax] = r1 * r2 * int(np.rint(nmax/r1/r2))
#        ms[icy(imax,-1)] = ms[imax]//r1
#        ms[icy(imax,1)] = ms[imax]//r2
#    print ms, 'ms'
#    return ms
#                      
#def getkpts_vasp(dir):
#    file1 = open(maindir+dir+'/'+kptsfile,'r')
#    kpointsfile = file1.readlines()
#    file1.close()
#    kpts_vasp = np.array(kpointsfile[3].split(),dtype=np.int16)
#    return kpts_vasp
#
#def writekpts_vasp(dir, mesh):
#    '''Write mesh m's to kpoints file, replacing previous mesh'''
#    file1 = open(maindir+dir+'/'+kptsfile,'r')
#    kpointsfile = file1.readlines()
#    file1.close
#    file2 = open(maindir+dir+'/'+kptsfile,'w')
#    kpointsfile[0] = 'Mesh generated by s/v and integer relations method  \n'    
#    kpointsfile[3] = '  '.join([str(mesh[i]) for i in [0,1,2]])+'\n'
#    file2.writelines(kpointsfile) 
#    file2.close()
#    return 


def writekpts_vasp_n(dir, mesh,n,type):
    '''Write mesh m's to kpoints file, using integer division for cubic and fcc meshes'''   
    file1 = open(maindir+dir+'/'+'KPOINTS','w')
    kpointsfile = []
    file1.close
    file2 = open(maindir+dir+'/'+kptsfile,'w')
    kpointsfile.append('For testing convergence vs n \n')
    kpointsfile.append('0 \n')   
    kpointsfile.append('Cartesian \n')
    if type == 'cubic':
        b = 1.0/n
        kpointsfile.append('%12.8f 0.0 0.0\n' % b)
        kpointsfile.append('0.0 %12.8f 0.0\n' % b)
        kpointsfile.append('0.0 0.0 %12.8f\n' %  b)
    if type == 'fcc':
        b = 1.0/n/2.0
        kpointsfile.append('0.0 %12.8f %12.8f\n' % (b,b))
        kpointsfile.append('%12.8f 0.0 %12.8f\n' % (b,b))
        kpointsfile.append('%12.8f %12.8f 0.0\n' % (b,b))
    else:
        sys.exit('Stop: type not found in writekpts_vasp_n')
    kpointsfile.append('0.5 0.5 0.5\n' %b) #shift
    file1.writelines(kpointsfile) 
    file1.close()
    return 

#
#def writejobfile(maindir,dir):
#    '''read from a template in maindir, and put dir in job name'''
#    file1 = open(maindir+'vaspjob','r')
#    jobfile = file1.readlines()
#    file1.close
#    for i in range(len(jobfile)):
#        jobfile[i]=jobfile[i].replace('myjob', dir)
#    file2 = open(maindir+dir+'/'+'vaspjob','w')
#    file2.writelines(jobfile) 
#    file2.close()
#    return 


################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50/AlIr/'
testfile = 'POSCAR.orig'
kptsfile = 'KPOINTS'
Nkppra = 10000

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
        print
        currdir = maindir + dir+'/'
        print dir + "************"
        file1 = open(maindir+dir+'/'+testfile,'r')
        poscar = file1.readlines()
        file1.close()
        if len(poscar) >0:
            os.chdir(currdir)
#            scale = np.sum(np.array(float(poscar[1])))
#            N = np.rint(Nkppra/np.sum(np.array(poscar[5].split(),dtype=np.int16))).astype(int) # number of kpts desired
#            reallatt[0,:] = np.array(poscar[2].split())
#            reallatt[1,:] = np.array(poscar[3].split())
#            reallatt[2,:] = np.array(poscar[4].split())
#            reallatt = scale*reallatt.astype(np.float)        
#            reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))

            writekpts_vasp_n(dir, mesh,n,type):
            writejobfile(maindir,dir)
            os.system('rm slurm*')
#            subprocess.call(['rm', 'slurm*'])
            subprocess.call(['rm', 'vasp.out'])
            subprocess.call(['rm', 'OUTCAR'])            
            subprocess.call(['cp','POSCAR.orig','POSCAR'])
#            subprocess.call(['sbatch', 'vaspjob'])
            os.chdir(maindir)
            
            # submit vasp job                      
print 'Done'


