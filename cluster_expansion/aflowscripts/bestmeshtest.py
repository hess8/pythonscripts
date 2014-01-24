#!/usr/bin/python
'''    Tests routine for finding best mesh via symmetry eigenvectors, for each structure in dir
'''
   
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64
from numpy.linalg import norm
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar
from bestmeshIter import bestmeshIter
fprec=float64

################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/AlIr/'
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr/'
testfile = 'POSCAR'
Nkppra = 10000

#reallatt = zeros((3,3))
os.chdir(maindir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
file1 = open('meshsummary.csv','w')
file1.write('Structure' + ',' +  'Lattice' + ',' + 'amax/amin' + ',' + 'Packing lattice' + ',' + 'Packing mesh' + ',' \
             + 'Improvement' + ',' + 'Mesh type' + ',' + 'Nmesh' + ',' + 'TargetNmesh' + '\n')
for dir in dirs:
    if testfile in os.listdir(dir):
        print
        print dir + '========================='
        path = maindir+dir+'/'
        os.chdir(path)
#        print readposcar('POSCAR',path)
        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = readposcar('POSCAR',path) #
        os.chdir(maindir)
#        print 'reciprocal lattice vectors (rows)';print reciplatt
        totatoms = sum(natoms)
        Nmesh = Nkppra/totatoms
        [meshvecs, Nmesh, targetNmesh, lattype, pfB, pfmax, meshtype] = bestmeshIter(reciplatt,Nmesh)
        pfimprove = round(pfmax/pfB , 1)
        a0 = norm(reciplatt[:,0]); a1 = norm(reciplatt[:,1]); a2 = norm(reciplatt[:,2]); 
        amax = max([a0,a1,a2]); amin =  min([a0,a1,a2])
        aratio = round(amax/amin, 1)
        file1.write(dir + ',' +  lattype + ',' + str(aratio) + ',' +str(pfB) + ',' + str(pfmax) + ',' + str(pfimprove) + ',  ' + meshtype + ','\
                     + str(Nmesh)+ ',' + str(targetNmesh)+  '\n')
file1.close()
        
print 'Done'