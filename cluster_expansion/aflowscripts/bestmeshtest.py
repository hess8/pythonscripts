#!/usr/bin/python
'''    Tests routine for finding best mesh via symmetry eigenvectors, for each structure in dir
'''
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint
from numpy.linalg import norm
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar
#from bestmeshIterJan22 import bestmeshIter
from bestmeshIter import bestmeshIter
fprec=float64

################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test101x/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test/'
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr/'
testfile = 'POSCAR'
Nkppra = 10000*10000

#reallatt = zeros((3,3))
os.chdir(maindir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
file1 = open('meshsummary.csv','w')
file1.write('Structure,Lattice,amax/amin,pfB,pf_orth,pf_orth2fcc,pf_maxpf, pf_pf2fcc, pfmax, meshtype' + ',' \
             + 'Improvement,fcc compatibility,Nmesh,TargetNmesh,Nmesh/Target,cbest' + '\n')
#for i,dir in enumerate(dirs):    
#    if testfile in os.listdir(dir) and i >8638:
for dir in dirs:
    if testfile in os.listdir(dir):        
        print 
        print dir + '=========================================================='
        path = maindir+dir+'/'
        os.chdir(path)
#        print readposcar('POSCAR',path)
        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = readposcar('POSCAR',path) #
        os.chdir(maindir)
#        print 'reciprocal lattice vectors (rows)';print reciplatt
        totatoms = sum(natoms)
        Nmesh = Nkppra/totatoms
        [meshvecs, Nmesh, targetNmesh, lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype,cbest, status] = bestmeshIter(reciplatt,Nmesh)
#        [K.vecs, K.Nmesh, B.Nmesh, B.lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype(B),status]
        pfimprove = round(pfmax/pfB , 1)
        a0 = norm(reciplatt[:,0]); a1 = norm(reciplatt[:,1]); a2 = norm(reciplatt[:,2]); 
        amax = max([a0,a1,a2]); amin =  min([a0,a1,a2])
        aratio = round(amax/amin, 1)
        format = 16*'%s,'+'%s'
        file1.write(format %  \
        (dir, lattype, str(aratio), str(pfB), str(pf_orth), str(pf_orth2fcc), str(pf_maxpf),str(pf_pf2fcc), str(pfmax), \
         meshtype, str(pfimprove), str(fcctype), str(rint(Nmesh)), str(targetNmesh), str(round(Nmesh/targetNmesh,3)),str(cbest),status+'\n'))
file1.close()
        
print 'Done'