#!/usr/bin/python
'''    Tests routine for finding best mesh via symmetry eigenvectors, for each structure in dir
'''
   
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64
import kmeshroutines as km
from kmeshroutines import nstrip
from bestmeshEigenIter import bestmeshEigenIter
fprec=float64

################# script #######################

#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr/'
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test/'
testfile = 'POSCAR'
Nkppra = 10000

#reallatt = zeros((3,3))
os.chdir(maindir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
for dir in dirs:
    if testfile in os.listdir(dir):
#        print
        print
        print dir + '========================='
        path = maindir+dir+'/'
        os.chdir(path)
#        print km.readposcar('POSCAR',path)
        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = km.readposcar('POSCAR',path) #
        os.chdir(maindir)
#        print 'reciprocal lattice vectors (rows)';print reciplatt
        totatoms = sum(natoms)
        Nmesh = Nkppra/totatoms
        bestmeshEigenIter(reciplatt,Nmesh)
print 'Done'