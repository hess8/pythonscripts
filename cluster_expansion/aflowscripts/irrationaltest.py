#!/usr/bin/python
''' tests whether the mesh numbers from the s/v method have irrational relationships:  sqrt 2, sqrt 3, sqrt 5 
Reads POSCAR info from aflow.in
'''
    
import sys,os
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
        print dir
        file1 = open(maindir+dir+'/'+testfile,'r')
        aflowin = file1.readlines()
        file1.close()
        for i,line in enumerate(aflowin):
            if 'VASP_POSCAR_MODE_EXPLICIT' in line:
                istart = i+1
                break #take only first instance (should be only one)
        cryststruc = np.array(aflowin[istart+2].split(), dtype=np.float)
#        print cryststruc
        reallatt =  km.lattice_vecs(cryststruc)
        reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
        natoms = np.array(aflowin[istart+3].split(),dtype=np.int16)
        print natoms
        totatoms=np.sum(natoms)
        positions = np.zeros((totatoms,3),dtype=np.float)
        postype = aflowin[istart+4].split()[0] #Direct or Cartesian
        where = 0
        for natom in natoms:
            for i in range(natom):
                for k in [0,1,2]:
                    positions[where,k] = float(aflowin[istart+5+where].split()[k])
                where += 1
        print positions
        totatoms=np.sum(natoms)
        km.create_poscar('From aflow.in, bch',1.0,reallatt,natoms,postype,positions)
        
        
        N = np.rint(Nkppra/totatoms).astype(int)      
#        [mesh_ns, irrat] = km.svmesh(N,reciplatt)
#        print mesh_ns, 's/v method'
#        if len(irrat)>0:
#            print dir, irrat




#                        
print 'Done'