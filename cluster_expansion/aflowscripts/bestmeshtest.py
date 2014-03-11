#!/usr/bin/python
'''    Tests routine for finding best mesh via symmetry eigenvectors, for each structure in dir
'''
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint
from numpy.linalg import norm
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar,create_poscar
#from bestmeshIterJan22 import bestmeshIter
from bestmeshIter import bestmeshIter
from bestmeshIter_vary_pf import bestmeshIter_vary_pf
fprec=float64

################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test101x/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/AlIr/'
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2/'
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2.10xNk/'

#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr34-50/'

testfile = 'POSCAR'
Nkppra = 10000*10

#reallatt = zeros((3,3))
os.chdir(maindir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
file1 = open('meshsummary.csv','a')
file1.write('Structure,Lattice,amax/amin,pfB,pf_orth,pf_orth2fcc,pf_maxpf, pf_pf2fcc, pfmax, meshtype' + ',' \
             + 'Improvement,fcc compatibility,Nmesh,TargetNmesh,Nmesh/Target,cbest' + '\n')
#for i,directory in enumerate(dirs):    

for directory in dirs:
    path = maindir+directory+'/'
    if testfile in os.listdir(path):        
        print 
        print directory + '=========================================================='
        os.chdir(path)
#        print readposcar('POSCAR',path)
        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = readposcar('POSCAR',path) #
        create_poscar('POSCAR',descriptor, scale, latticevecs, natoms, postype, positions, path) #just to remove the scale problem
        os.chdir(maindir)
#        print 'reciprocal lattice vectors (rows)';print reciplatt
        totatoms = sum(natoms)
        Nmesh = Nkppra/totatoms
#RESTORE THIS        [meshvecs, Nmesh, targetNmesh, lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype,cbest, status] = bestmeshIter(reciplatt,Nmesh)
 
# for trials where we want to vary pf for testing       
        meshesfile = open('meshesfile','w')
        meshesfile.write(directory+' ============\n')
        meshesfile.close()
        [meshvecs, Nmesh, targetNmesh, lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype,cbest, status] = bestmeshIter_vary_pf(reciplatt,Nmesh,path)
# End trials

#        [K.vecs, K.Nmesh, B.Nmesh, B.lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype(B),status]
        pfimprove = round(pfmax/pfB , 1)
        a0 = norm(reciplatt[:,0]); a1 = norm(reciplatt[:,1]); a2 = norm(reciplatt[:,2]); 
        amax = max([a0,a1,a2]); amin =  min([a0,a1,a2])
        aratio = round(amax/amin, 1)
        format = 16*'%s,'+'%s'
        file1.write(format %  \
        (directory, lattype, str(aratio), str(pfB), str(pf_orth), str(pf_orth2fcc), str(pf_maxpf),str(pf_pf2fcc), str(pfmax), \
         meshtype, str(pfimprove), str(fcctype), str(rint(Nmesh)), str(targetNmesh), str(round(Nmesh/targetNmesh,3)),str(cbest),status+'\n'))
file1.close()
        
print 'Done'