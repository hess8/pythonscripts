#!/usr/bin/python
''' tests. '''
    
import sys,os
import numpy as np
################# functions #######################
def svmesh(N,vecs):
    '''N: points desired.  vecs the lattice vectors as numpy array (reciprocal in our thinking)
    output:  n1, n2, n3, the number of divisions along each RLV for the mesh'''
    u = np.linalg.norm(np.cross(vecs[0,:],vecs[1,:]))
    v = np.linalg.norm(np.cross(vecs[1,:],vecs[2,:]))
    w = np.linalg.norm(np.cross(vecs[2,:],vecs[0,:]))
    n1 = N**(1/3.0) * u**(1/3.0) * w**(1/3.0) / v**(2/3.0)
    n2 = N**(1/3.0) * u**(1/3.0) * v**(1/3.0) / w**(2/3.0)
    n3 = N**(1/3.0) * v**(1/3.0) * w**(1/3.0) / u**(2/3.0)
    print 'order type', ordertype([0,n1,n2,n3])
#    print n1,n2,n3, n1*n2*n3

    p = n2/n1
    q = n3/n2
    r = n1/n3
    delta = 10**-6
    PQR = np.array([0,0,0]).astype(int)
    '''   Define triangle as 
               m1
            P      R
        
        m2     Q     m3
    The integer relations (P,Q,R) 'point' to the larger integer.  If they are CCW, record +1.  
    If CW, -1 
                  
           ''' 
   
    if abs(np.rint(p)-p)<delta:
        PQR[0] = 1
    elif abs(np.rint(1/p)-(1/p))<delta:
        PQR[0] = -1
    if abs(np.rint(q)-q)<delta:
        PQR[1] = 1
    elif abs(np.rint(1/q)-(1/q))<delta:
        PQR[1] = -1    
    if abs(np.rint(r)-r)<delta:
        PQR[2] = 1
    elif abs(np.rint(1/r)-(1/r))<delta:
        PQR[2] = -1
   
    Nrels = np.sum(np.abs(PQR)) #number of integer relations)
    print 'Nrels', Nrels, PQR          
    
    return np.rint(np.array([n1,n2,n3])).astype(int) #round then convert to integers

def getkpts_vasp(dir):
    file1 = open(maindir+dir+'/'+kptsfile,'r')
    kpointsfile = file1.readlines()
    file1.close()
    kpts_vasp = np.array(kpointsfile[3].split(),dtype=np.int16)
    return kpts_vasp

def ordertype(nlist):
    '''defining 6 orderings, so can reverse later.  I'm sure there is a better way to do this'''
    if  (nlist[1] >= nlist[2]) & (nlist[2] >= nlist[3]):
        return 'a'
    if  (nlist[1] >= nlist[3]) & (nlist[3] >= nlist[2]):
        return 'b'
    if  (nlist[2] >= nlist[1]) & (nlist[1] >= nlist[3]):
        return 'c'
    if  (nlist[2] >= nlist[3]) & (nlist[3] >= nlist[1]):
        return 'd'
    if  (nlist[3] >= nlist[1]) & (nlist[1] >= nlist[2]):
        return 'e'
    if  (nlist[3] >= nlist[2]) & (nlist[2] >= nlist[1]):
        return 'f'
  
   
################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAktest/AlIr/'
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
            print 'natoms', np.array(poscar[5].split(),dtype=np.int16)
            N = np.rint(Nkppra/np.sum(np.array(poscar[5].split(),dtype=np.int16))).astype(int) # number of kpts desired
            reallatt[0,:] = np.array(poscar[2].split())
            reallatt[1,:] = np.array(poscar[3].split())
            reallatt[2,:] = np.array(poscar[4].split())
            reallatt = scale*reallatt.astype(np.float)        
#            print reallatt
            reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
            print reciplatt
            mesh_ns = svmesh(N,reciplatt)
            print mesh_ns
            print getkpts_vasp(dir), 'Aflow'
                        
print 'Done'