#!/usr/bin/python
''' tests whether the mesh numbers from the s/v method have irrational relationships:  sqrt 2, sqrt 3, sqrt 5 
Reads POSCAR info from aflow.in
'''
    
import sys,os,subprocess
import numpy as np
import kmeshroutines as km
from kmeshroutines import nstrip
from bestmesh import bestmesh
#from poscar import POSCAR as POSCAR 

def check_out(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    output = process.communicate()
    retcode = process.poll()
    if retcode:
            raise subprocess.CalledProcessError(retcode, command, output=output[0])
    return output          
   
################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr/'
testfile = 'POSCAR'
Nkppra = 10000

#reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
#        print
        print
        print dir + '========================='
        path = maindir+dir+'/'
        os.chdir(path)
        print km.readposcar('POSCAR',path)
        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = km.readposcar('POSCAR',path) #
        os.chdir(maindir) 
        totatoms = np.sum(natoms)
        Nmesh = Nkppra/totatoms
        bestmesh(reciplatt,Nmesh)

       

#        print 'lattice from aconvasp --sprim < POSCAR0 > POSCAR'
##        print 'lattice from aconvasp --sprim < POSCAR0 > POSCAR'
#        print latticevecs
#        print
#        print 'reciprocal lattice vectors'
#        print reciplatt
#        print
#        reciplatt = np.array((pos.bvecs[0],pos.bvecs[1],pos.bvecs[2]))
#        print reciplatt
##        print 'bvecs'
##        print pos.bvecs
##        print pos.bvecs[1,:]   
#        [mesh_ns, irrat] = km.svmesh(N,pos.bvecs)        
#        [mesh_ns, irrat] = km.svmesh(N,reciplatt)
##        km.writekpts_vasp(dir,mesh_ns) #correct kmesh
#        print mesh_ns, 's/v method'
#        if len(irrat)>0:
#            print dir, irrat 
                        
 
print 'Done'