#!/usr/bin/python
'''    Tests routine for finding best mesh via symmetry eigenvectors, for each structure in dir
'''


import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy
fprec=float64

################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test101x/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2/f5411/'
maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2.10xNk/f5411/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test/f3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr34-50/'

testfile = 'POSCAR'

#reallatt = zeros((3,3))
os.chdir(maindir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
#file1 = open('varypf.csv','a')
#file1.write('Structure,Lattice,amax/amin,pfB,pf_orth,pf_orth2fcc,pf_maxpf, pf_pf2fcc, pfmax, meshtype' + ',' \
#             + 'Improvement,fcc compatibility,Nmesh,TargetNmesh,Nmesh/Target,cbest' + '\n')
#for i,directory in enumerate(dirs):    
print dirs
writeEnergiesOszicar(dirs) 
writedirnames(dirs)
writeNk(dirs)


################# summary #################
outfile = open('vary_pf.csv','w')
outfile.write('pf,energy\n')
 

#os.chdir(mainDir)
file = open(maindir+'names','r')
names = nstrip(file.readlines())
file.close()

file = open(maindir+'energies','r')
energies = nstrip(file.readlines())
file.close()

for i in range(len(names)):
    linei = names[i]+','+energies[i]+'\n'        
    outfile.write(linei)    
outfile.close()  

struct = maindir.split('/')[-2]
print struct

plotxy(names,energies,'vary_pf','Struct + Vasp energy vs packing fraction','Packing fraction','eV')

      
print 'Done'

