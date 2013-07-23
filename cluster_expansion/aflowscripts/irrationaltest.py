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
                print 'istart', istart
                break #take only first instance (should be only one)
           
            
        scale = np.sum(np.array(float(aflowin[istart+1])))
#            print 'natoms', np.array(poscar[5].split(),dtype=np.int16)
        N = np.rint(Nkppra/np.sum(np.array(aflowin[istart+5].split(),dtype=np.int16))).astype(int) # number of kpts desired
        reallatt[0,:] = np.array(aflowin[istart+2].split())
        reallatt[1,:] = np.array(aflowin[istart+3].split())
        reallatt[2,:] = np.array(aflowin[istart+4].split())
        reallatt = scale*reallatt.astype(np.float)        
#            print reallatt
        reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
#            print reciplatt
        mesh_ns = km.svmesh(N,reciplatt)
        print mesh_ns, 's/v method'


#                        
print 'Done'