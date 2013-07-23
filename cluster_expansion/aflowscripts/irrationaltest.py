#!/usr/bin/python
''' tests whether the mesh numbers from the s/v method have irrational relationships:  sqrt 2, sqrt 3, sqrt 5 
Reads POSCAR info from aflow.in
'''
    
import sys,os
import numpy as np
import kmeshroutines as km                
   
################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/AlIr/'
testfile = 'aflow.in'
Nkppra = 10000

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
        print
        print dir
        file1 = open(maindir+dir+'/'+testfile,'r')
        aflowin = file1.readlines()
        file1.close()
        for i,line in enumerate(aflowin):
            if 'VASP_POSCAR_MODE_EXPLICIT' in line:
                istart = i+1
                break #take only first instance (should be only one)
        cryststruc = np.array(aflowin[istart+2].split(), dtype=np.float)
        print cryststruc
        print km.lattice_vecs(cryststruc)
#        mesh_ns = km.svmesh(N,reciplatt)
#        print mesh_ns, 's/v method'


#                        
print 'Done'