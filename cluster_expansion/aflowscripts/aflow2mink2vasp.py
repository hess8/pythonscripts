#!/usr/bin/python
''' Reads POSCAR from aflow run.  Replaces lattice and atom positions with a lattice that has a
 mink reduced reciprocal lattice.  Overwrites the kpoints file.  
'''
    
import sys,os,subprocess
import numpy as np
import kmeshroutines as km
from kmeshroutines import nstrip
from poscar import POSCAR 
import poscar

def check_out(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    output = process.communicate()
    retcode = process.poll()
    if retcode:
            raise subprocess.CalledProcessError(retcode, command, output=output[0])
    return output          
   
################# script #######################

maindir = '/fslhome/bch/cluster_expansion/alir/testf10/AlIr/'
aflowdir = '/fslhome/bch/cluster_expansion/alir/testf10mirror/AlIr/'
testfile = 'aflow.in'
Nkppra = 10000

reallatt = np.zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
#        print
        print
        print dir + '========================='
        mirdir = aflowdir+dir+'/'
        path = maindir+dir+'/'
        os.chdir(path)
#        os.system('rm POSCAR*')
        os.system('cp %sKPOINTS .' % mirdir)
        os.system('cp %sINCAR .' % mirdir)
        os.system('cp %sPOTCAR .' % mirdir)
        os.system('cp %svaspjob .' % mirdir)
        os.system('cp %sPOSCAR POSCARaflow' % mirdir)        
#        totatoms = km.aflow2poscar(path)
#        os.system('aconvasp --sprim < POSCAR0 > POSCAR')
#        back = subprocess.check_output(["echo", "Hello World!"])  
#        back = check_out(["echo", "Hello World!"])
#        back = os.system('mink_reduction.py < POSCAR0 ')
#        back = subprocess.check_output(['mink_reduction.py',' <',' POSCAR0 '], shell=True,)
#        back = nstrip(subprocess.check_output(['mink_reduction.py < POSCAR0'], shell=True,).split())
#        print 'back'
#        print back

        pfile = open('POSCARaflow','r')
        rlines  = [i.strip() for i in pfile.readlines()]       
        pos = POSCAR(lines=rlines)
        totatoms = np.sum(pos.types)
        N = np.rint(Nkppra/totatoms).astype(int)   
        oldROD = pos.rod
        oldbvecs = np.array(pos.bvecs)
        reallatt = np.array((pos.avecs[0],pos.avecs[1],pos.avecs[2]))
        print "Real lattice"
        print  reallatt        
        print "Real lattice orth defect:",pos.orthogonality_defect
        print "Recp lattice orth defect:",pos.rod
        pos.kmink_reduce(1e-1) ####THIS DOES NOT CHANGE THE BVECS.  HAVE TO WRITE A NEW POSCAR, AS THE PROPERTIES ARE READ ONLY
        sys.stdout.flush()
        newROD = pos.rod
        newbvecs = np.array(pos.bvecs)
        reciplatt = np.array((pos.bvecs[0],pos.bvecs[1],pos.bvecs[2]))
#        print reciplatt
        test = []
        if abs(newROD - oldROD) < 1e-3 and np.equal(newbvecs,oldbvecs).all:
#        if abs(newROD - oldROD) < 1e-3 :
            print "Recip lattice already reduced"
            print reciplatt
            pos.name = pos.name + ' k-space mink reduced OK'    
        else:
            pos.name = pos.name + ' New lattice for k-space mink reduced'
            print "Old recip lattice"
            print oldbvecs
            print "New recip lattice orth defect", newROD
            print "New recip lattice"
            print newbvecs
            print "New real lattice"
            reallatt  = 2*np.pi**np.linalg.inv(np.transpose(reciplatt)) 
            print  reallatt
        pos.write_poscar('POSCAR')
#        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = km.readposcar('POSCAR0',path) #
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
        [mesh_ns, irrat] = km.svmesh(N,reciplatt)
        km.writekpts_vasp(maindir,dir,'KPOINTS',mesh_ns) #correct kmesh
        
        print mesh_ns, 's/v method'
#        if len(irrat)>0:
#            print dir, irrat
        print irrat
        os.system('rm slurm*')
        subprocess.call(['rm', 'vasp.out'])
        subprocess.call(['rm', 'OUTCAR'])            
#        subprocess.call(['sbatch', 'vaspjob'])                        
        os.chdir(maindir)   
print 'Done'