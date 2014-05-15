#!/usr/bin/python
''' tests. '''
    
import sys,os
import numpy as np
import kmeshroutines as km
################# functions #######################
                   

#def writekpts_vasp(dir, mesh):
#    '''Write mesh m's to kpoints file, replacing previous mesh'''
#    file1 = open(maindir+dir+'/'+kptsfile,'r')
#    kpointsfile = file1.readlines()
#    file1.close
#    file2 = open(maindir+dir+'/'+kptsfile,'w')
#    kpointsfile[3] = '  '.join([str(mesh[i]) for i in [0,1,2]])+'\n'
#    file2.writelines(kpointsfile) 
#    file2.close()
#    return 

   
################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50/AlIr/'
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
#            print 'natoms', np.array(poscar[5].split(),dtype=np.int16)
            N = np.rint(Nkppra/np.sum(np.array(poscar[5].split(),dtype=np.int16))).astype(int) # number of kpts desired
            reallatt[0,:] = np.array(poscar[2].split())
            reallatt[1,:] = np.array(poscar[3].split())
            reallatt[2,:] = np.array(poscar[4].split())
            reallatt = scale*reallatt.astype(np.float)        
#            print reallatt
            reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
#            print reciplatt
            mesh_ns = km.svmesh(N,reciplatt)
            print mesh_ns, 's/v method'
            if 'KPOINTS.relax1' in os.listdir(dir):
                os.chdir(dir)
                os.system('cp KPOINTS.relax1 KPOINTS')
                os.chdir(maindir)               
            print km.getkpts_vasp(dir), 'Aflow'
#            writekpts_vasp(dir, mesh_ns)
#                        
print 'Done'