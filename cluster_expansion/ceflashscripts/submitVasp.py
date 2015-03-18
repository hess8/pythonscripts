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

################# script #######################

#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50/AlIr/'
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2/f5179/'
testfile = 'POSCAR'

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
        currdir = maindir + dir+'/'
        print dir #+ "************"
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

            os.system('rm slurm*')
#            subprocess.call(['rm', 'slurm*'])
            subprocess.call(['rm', 'CHG*'])
            subprocess.call(['rm', 'OUTCAR'])            

            subprocess.call(['sbatch', 'vaspjob'])
            os.chdir(maindir)
            
            # submit vasp job                      
print 'Done'



