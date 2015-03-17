#!/usr/bin/python
''' tests whether the mesh numbers from the s/v method have irrational relationships:  sqrt 2, sqrt 3, sqrt 5 
Reads POSCAR info from aflow.in
'''
    
import sys,os,subprocess
import numpy as np
import kmeshroutines as km
from kmeshroutines import nstrip
from poscar import POSCAR as POSCAR 

def check_out(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    output = process.communicate()
    retcode = process.poll()
    if retcode:
            raise subprocess.CalledProcessError(retcode, command, output=output[0])
    return output          
   
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
        print
        print dir + '========================='
        path = maindir+dir+'/'
        os.chdir(path)
#        os.system('rm POSCAR*')
        os.system('cp ../../KPOINTS .')
        os.system('cp ../../INCAR .')
        os.system('cp ../../POTCAR .')        
        totatoms = km.aflow2poscar(path)
#        os.system('aconvasp --sprim < POSCAR0 > POSCAR')
#        back = subprocess.check_output(["echo", "Hello World!"])  
#        back = check_out(["echo", "Hello World!"])
#        back = os.system('mink_reduction.py < POSCAR0 ')
#        back = subprocess.check_output(['mink_reduction.py',' <',' POSCAR0 '], shell=True,)
#        back = nstrip(subprocess.check_output(['mink_reduction.py < POSCAR0'], shell=True,).split())
#        print 'back'
#        print back
        N = np.rint(Nkppra/totatoms).astype(int)
        pfile = open('POSCAR0','r')
        rlines  = [i.strip() for i in pfile.readlines()]       
        pos = POSCAR(lines=rlines)   
        oldROD = pos.rod
        oldbvecs = np.array(pos.bvecs)
        print "Real lattice orth defect:",pos.orthogonality_defect
        print "Recp lattice orth defect:",pos.rod
        pos.kmink_reduce(1e-1)
        sys.stdout.flush()
        newROD = pos.rod
        newbvecs = np.array(pos.bvecs)
        reciplatt = np.array((pos.bvecs[0],pos.bvecs[1],pos.bvecs[2]))
#        print reciplatt
        test = []
        if abs(newROD - oldROD) < 1e-3 and np.equal(newbvecs,oldbvecs).all:
#        if abs(newROD - oldROD) < 1e-3 :
            print "Recip lattice already reduced"         
        else:
            print "Old recip lattice"
            print oldbvecs
            print "New recip lattice orth defect", newROD
            print "New recip lattice"
            print newbvecs
            print "New real lattice"
            reallatt  = 2*np.pi**np.linalg.inv(np.transpose(reciplatt)) 
            print  reallatt

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
#        km.writekpts_vasp(dir,mesh_ns) #correct kmesh
        print mesh_ns, 's/v method'
        if len(irrat)>0:
            print dir, irrat                        
        os.chdir(maindir)   
print 'Done'