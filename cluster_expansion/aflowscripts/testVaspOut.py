#/usr/bin/env python
import time, os, subprocess, sys
'''goes through every subdir and checks files for errors using grep, etc
'''
import numpy as np   

import kmeshroutines as km


################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50corr/AlIr/'
testfile = 'vasp.out'
kptsfile = 'KPOINTS'

os.chdir(maindir)
print os.getcwd()
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
#standardout = sys.__stdout__
#sys.stdout = open('testvasp.dat', 'w')
for dir in dirs:
    print
    print dir +'   ***************'
    if testfile in os.listdir(dir):
        currdir = maindir + dir +'/'
        os.chdir(currdir)
        print currdir
        file1 = open(testfile,'r')
        vaspout = file1.readlines()
        file1.close()
        for i in [4,5,6,7]:
            print vaspout[i].replace('\n','')
        kmesh = km.getkpts_vasp(currdir)
        print kmesh, kmesh[0]*kmesh[1]*kmesh[2] 
        [natoms,reallatt,reciplatt] = km.readposcar(currdir) 
        Nkppra = 10000
        N = np.rint(Nkppra/natoms).astype(int)
        print 'natoms', natoms
        print 'reciplatt'
        print reciplatt
        print 'mesh', km.svmesh(N,reciplatt)
#        os.system('grep -i "bad news" slurm*')
#        os.system('grep -i kpts slurm*')
#        os.system('grep -i exceed slurm*')
#        os.system('grep -i sgrcon slurm*') 
#        os.system('grep -i bad vasp.out')
#        print subprocess.call(['grep','-i','bad','vasp.out']) 
#        slurm = os.system('find slurm*')
        slurms = [f for f in os.listdir(os.getcwd()) if 'slurm' in f]
        for slurm in slurms:
            km.regpy_nocase('bad',currdir+slurm)
#        print subprocess.check_output(['grep','-i','bad','vasp.out'])
#        print slurm
#        km.regpy_nocase('bad',currdir+ slurm)
        km.regpy_nocase('bad',currdir+'vasp.out')
#        subprocess.call(['grep','-i','bad','vasp.out']) 
#        sys.stdout.flush()                   
        os.chdir(maindir)   
#sys.__stdout__ =    standardout                     
print 'Done'



