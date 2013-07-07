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
    print n1,n2,n3, n1*n2*n3
    return np.rint(np.array([n1,n2,n3])).astype(int) #round then convert to integers
    
################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA4Gb4proc2/AlIr/'
testfile = 'POSCAR'
N =100

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs=[d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
for dir in dirs:
    if testfile in os.listdir(dir):
        print dir
        file1 = open(maindir+dir+'/'+testfile,'r')
        poscar = file1.readlines()
        if len(poscar) >0:
            scale = float(poscar[1])
            reallatt[0,:] = np.array(poscar[2].split())
            reallatt[1,:] = np.array(poscar[3].split())
            reallatt[2,:] = np.array(poscar[4].split())
            reallatt = scale*reallatt.astype(np.float)         
#            print reallatt
            reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
            print reciplatt
            mesh_ns = svmesh(N,reciplatt)
            print mesh_ns
        
print 'Done'