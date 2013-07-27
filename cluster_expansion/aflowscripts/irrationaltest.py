#!/usr/bin/python
''' tests whether the mesh numbers from the s/v method have irrational relationships:  sqrt 2, sqrt 3, sqrt 5 
Reads POSCAR info from aflow.in
'''
    
import sys,os,subprocess
import numpy as np
import kmeshroutines as km                
   
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
#        print dir
        path = maindir+dir+'/'
        totatoms = km.aflow2poscar(path)
        os.chdir(path)
        os.system('aconvasp --sprim < POSCAR0 > POSCAR')
#        subprocess.call(['aconvasp','--sprim','<','POSCAR0','>','POSCAR'])
        os.chdir(maindir)       
        N = np.rint(Nkppra/totatoms).astype(int)
        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = km.readposcar('POSCAR',path) #
#        print [descriptor, scale, latticevecs, natoms]
        [mesh_ns, irrat] = km.svmesh(N,reciplatt)
#        print mesh_ns, 's/v method'
        if len(irrat)>0:
            print dir, irrat
#                        
print 'Done'