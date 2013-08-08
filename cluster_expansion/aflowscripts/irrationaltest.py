#!/usr/bin/python
''' tests whether the mesh numbers from the s/v method have irrational relationships:  sqrt 2, sqrt 3, sqrt 5 
Reads POSCAR info from aflow.in
'''
    
import sys,os,subprocess
import numpy as np
import kmeshroutines as km
from poscar import POSCAR as POSCAR          
   
################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50b/AlIr/'
testfile = 'aflow.in'
Nkppra = 10000

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
#        print
        print
        print dir + '========================='
        path = maindir+dir+'/'
        os.chdir(path)
#        os.system('rm POSCAR*')
        totatoms = km.aflow2poscar(path)
#        os.system('aconvasp --sprim < POSCAR0 > POSCAR')
    
        N = np.rint(Nkppra/totatoms).astype(int)
        pfile = open('POSCAR0','r')
        rlines  = [i.strip() for i in pfile.readlines()]       
        pos = POSCAR(lines=rlines)       
        print "real  od:",pos.orthogonality_defect
        print "recip od:",pos.rod
#        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = km.readposcar('POSCAR0',path) #
#        print 'lattice from aconvasp --sprim < POSCAR0 > POSCAR'
##        print 'lattice from aconvasp --sprim < POSCAR0 > POSCAR'
#        print latticevecs
#        print
#        print 'reciprocal lattice vectors'
#        print reciplatt
        print
        reciplatt = np.array((pos.bvecs[0],pos.bvecs[1],pos.bvecs[2]))
        print reciplatt
#        print 'bvecs'
#        print pos.bvecs
#        print pos.bvecs[1,:]   
#        [mesh_ns, irrat] = km.svmesh(N,pos.bvecs)        
        [mesh_ns, irrat] = km.svmesh(N,reciplatt)
        print mesh_ns, 's/v method'
        os.chdir(maindir)   
        if len(irrat)>0:
            print dir, irrat                        
print 'Done'