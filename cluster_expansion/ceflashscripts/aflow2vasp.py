#!/usr/bin/python
''' Assumes vasp is ready to do mink reduction and chooses mesh. This overwrites puts "Minkowski Monkhorst-Pack" in place of "Monkhorst-Pack" 
in KPOINTS.  Replaces nx, ny, nz with "number of kpts per reciprocal atom"
'''
import sys,os,subprocess,time
import numpy as np
from numpy import pi
import kmeshroutines as km
from kmeshroutines import nstrip
#from poscar import POSCAR 
#import poscar

def check_out(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    output = process.communicate()
    retcode = process.poll()
    if retcode:
            raise subprocess.CalledProcessError(retcode, command, output=output[0])
    return output          
  
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50d.clean/AlIr/';vaspexec = 'vasp533' 
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50d/AlIr/';vaspexec = 'vasp533mod'
testfile = 'aflow.in'
Nkppra = 10000

reallatt = np.zeros((3,3))

os.chdir(maindir)
os.system('rm ../slurm*')
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
#        print
        print
        print dir + '========================='
#        mirdir = aflowdir+dir+'/'
        path = maindir+dir+'/'
        print path
        os.chdir(path)
#        os.system('cp %sKPOINTS .' % maindir)
#        os.system('cp %sINCAR .' % maindir)
#        os.system('cp %sPOTCAR .' % maindir)
        os.system('cp ../../KPOINTS .') #standard files 2 levels up
        os.system('cp ../../INCAR .')
        os.system('cp ../../POTCAR .')
        if vaspexec == 'vasp533mod':
            km.writekpts_vasp(maindir,dir+'/','KPOINTS',Nkppra) #correct 2 lines
            nameadd = 'mod'
        elif vaspexec == 'vasp533':
            nameadd = 'clean'
        km.writejobfile(maindir,dir+'/','vaspjob',nameadd,vaspexec)
#        if len(irrat)>0:
#            print dir, irrat
#        print irrat
        os.system('rm slurm*')
#        os.system('cp POSCAR.orig POSCAR') #in case aflow setup went too far. 
        km.aflow2poscar(path) #generates POSCAR from aflow.in.  A
        subprocess.call(['rm', 'vasp.out'])
        subprocess.call(['rm', 'OUTCAR'])            
        subprocess.call(['sbatch', 'vaspjob'])                        
        os.chdir(maindir)   
#km.checkq('bch') #loops for days
print 'Done'