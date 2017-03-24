#!/usr/bin/python
import time, os, subprocess
'''For each dir in maindir: copies vaspinput, creates kpoints file with correct mesh 
(just relative to the recip lattice),reads a jobfiles from the maindir,
writes the structure number etc to the job name, and submits a vasp job.  
Assumes we already have POSCAR from aconvaspPoscar.py (and before that train_structs2aflow.py)
'''
    
import sys,os
import numpy as np
#from kmeshroutines.py import *



def writekpts_vasp_n(path,n,type):
    '''Write mesh vectors to kpoints file, using integer division for cubic and fcc meshes'''   
    file1 = open(path +'KPOINTS','w')
    kpointsfile = []
    if type == 'cubic': kpointsfile.append('%i kpoints for cubic integer n=%i\n' %(n**3,n))
    if type == 'fcc': kpointsfile.append('%i kpoints for fcc integer n=%i\n' %(4*n**3,n))
    kpointsfile.append('0 \n')   
    kpointsfile.append('Reciprocal \n')
    print type
    if type == 'cubic':
        b = 1.0/n
        kpointsfile.append('%12.8f 0.0 0.0\n' % b)
        kpointsfile.append('0.0 %12.8f 0.0\n' % b)
        kpointsfile.append('0.0 0.0 %12.8f\n' %  b)
    elif type == 'fcc':
        b = 1.0/n/2.0
        kpointsfile.append('0.0 %12.8f %12.8f\n' % (b,b))
        kpointsfile.append('%12.8f 0.0 %12.8f\n' % (b,b))
        kpointsfile.append('%12.8f %12.8f 0.0\n' % (b,b))
    else:
        sys.exit('Stop: type not found in writekpts_vasp_n')
#     kpointsfile.append('0.5 0.5 0.5\n' ) #shift
    kpointsfile.append('0.0 0.0 0.0\n' ) #shift
    file1.writelines(kpointsfile) 
    file1.close()
    return 
 
def writejobfile(path,n,type):
    '''read from a template in maindir, and put dir in job name'''
    file1 = open(path +'vaspjob','r')
    jobfile = file1.readlines()
    file1.close
    for i in range(len(jobfile)):
        jobfile[i]=jobfile[i].replace('myjob', dir+'_%i_%s' %(n,type))
    file2 = open(path+'/'+'vaspjob','w')
    file2.writelines(jobfile) 
    file2.close()
    return 

def createdir(path,n,type):
    newdir = path + '%s_%i/' % (type,n)
    if not os.path.isdir(newdir):
        os.system('mkdir %s' % newdir)
#    os.chdir(newdir)
    os.system ('cp %s* %s' % (vaspinputdir,newdir))
    os.system ('cp %sPOSCAR %s' % (path,newdir))  
    writekpts_vasp_n(newdir,n,type)
    writejobfile(newdir,n,type)  
    return


################# script #######################

maindir = '/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/cubicTest/'
type = 'cubic'
testfile = 'POSCAR'
vaspinputdir = '/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/vaspinputShort/'
# Nkppra = 10000

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
        print
        currdir = maindir + dir+'/'
        print dir + "************"
        file1 = open(currdir+testfile,'r')
        poscar = file1.readlines()
        file1.close()
        if len(poscar) > 0:
            os.chdir(currdir)
#            scale = np.sum(np.array(float(poscar[1])))
#            N = np.rint(Nkppra/np.sum(np.array(poscar[5].split(),dtype=np.int16))).astype(int) # number of kpts desired
#            reallatt[0,:] = np.array(poscar[2].split())
#            reallatt[1,:] = np.array(poscar[3].split())
#            reallatt[2,:] = np.array(poscar[4].split())
#            reallatt = scale*reallatt.astype(np.float)        
#            reciplatt = 2*np.pi*np.transpose(np.linalg.inv(reallatt))
          
            

            os.system('rm slurm*')
            subprocess.call(['rm', 'vasp.out'])
            subprocess.call(['rm', 'OUTCAR'])          
#            subprocess.call(['cp','POSCAR.orig','POSCAR'])
#            subprocess.call(['sbatch', 'vaspjob'])
            
            # Now create new dirs with different characteristics
            type = 'cubic'
            for n in range(1,24,1):
                createdir(currdir,n,type)
#             type = 'fcc'
#             for n in range(1,15):
#                 createdir(currdir,n,type)  
#            type = 'bcc'
#            for n in range(1,18):
#                createdir(currdir,n,type)                       
            # submit vasp job 
        newdirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]) 
        for newdir in newdirs:
            os.chdir(newdir)
            subprocess.call(['sbatch', 'vaspjob'])
            os.chdir(currdir)
        os.chdir(maindir)                 
print 'Done'



