#!/usr/bin/env python
'''    Tests routine for finding best mesh for each structure in dir
'''
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint
from numpy.linalg import norm
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar,create_poscar
#from bestmeshIterJan22 import bestmeshIter
from bestmeshIter import bestmeshIter
from bestmeshIter_vary_pf import bestmeshIter_vary_pf
from bestmeshIter_vary_N import bestmeshIter_vary_N
fprec=float64
import meshConstruct

################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test101x/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test/'

#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2.10xNk/'

#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test.10xNk/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test.noshift/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test10^3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/f3varyN/'
# maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test/'
maindir = '/fslhome/bch/cluster_expansion/meshConstruct/AlAl'


#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/testSi/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/testMP/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr34-50/'
#maindir = '/fslhome/bch/cluster_expansion/sisi/test10^3/'
#maindir = '/fslhome/bch/cluster_expansion/sisi/test10^4/'

meshc = meshConstruct.meshConstruct() #instance

testfile = 'POSCAR'
# Nkppra = 10000#*10  
Nkppra = 1000#*10 
#reallatt = zeros((3,3))
os.chdir(maindir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
if os.path.exists('meshsummary.csv'):
    file1 = open('meshsummary.csv','a')
else: 
    file1 = open('meshsummary.csv','w')
file1.write('Structure,Lattice,amax/amin,pfB,pf_orth,pf_orth2fcc,pf_maxpf, pf_pf2fcc, pfmax, meshtype' + ',' \
             + 'Improvement,fcc compatibility,Nmesh,TargetNmesh,Nmesh/Target,cbest' + '\n')
#for i,dir in enumerate(dirs):    

for dir in dirs:
    path = maindir+'/'+dir
    print 'os.listdir(path)',os.listdir(path)
    if testfile in os.listdir(path):        
        print 
        print dir + '=========================================================='
        os.chdir(path)
#        print readposcar('POSCAR',path)
        [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = readposcar('POSCAR',path) #
        create_poscar('POSCAR',descriptor, scale, latticevecs, natoms, postype, positions, path) #just to remove the scale problem
        os.chdir(maindir)
#        print 'reciprocal lattice vectors (rows)';print reciplatt
        totatoms = sum(natoms)
        targetNmesh = Nkppra/totatoms
#RESTORE THIS        [meshvecs, Nmesh, targetNmesh, lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype,cbest, status] = bestmeshIter(reciplatt,Nmesh)
 
# for trials where we want to vary pf for testing       
#         meshesfile = open('meshesfile','w')
#         meshesfile.write(dir+' ============\n')
        [meshvecs, Nmesh, lattype, pfB, pf, status] = meshc.relaxMeshSym(reciplatt,targetNmesh,path)
#        bestmeshIter_vary_N(reciplatt,Nmesh,path)
# End trials

#        [K.vecs, K.Nmesh, B.Nmesh, B.lattype, pfB, pf_orth, pf_orth2fcc, pf_maxpf, pf_pf2fcc, pfmax, meshtype, fcctype(B),status]
        pfimprove = round(pfmax/pfB , 1)
        a0 = norm(reciplatt[:,0]); a1 = norm(reciplatt[:,1]); a2 = norm(reciplatt[:,2]); 
        amax = max([a0,a1,a2]); amin =  min([a0,a1,a2])
        aratio = round(amax/amin, 1)
        format = 16*'%s,'+'%s'
        file1.write(format %  \
        (dir, lattype, str(aratio), str(pfB), str(pf_orth), str(pf_orth2fcc), str(pf_maxpf),str(pf_pf2fcc), str(pfmax), \
         meshtype, str(pfimprove), str(fcctype), str(rint(Nmesh)), str(targetNmesh), str(round(Nmesh/targetNmesh,3)),str(cbest),status+'\n'))
        file2 = open('lattype','w'); file2.write(lattype); file2.close()
file1.close()
        
print 'Done'