
#!/usr/bin/python
'''    
'''

import sys,os,subprocess
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from numpy import zeros,transpose,array,sum,float64,rint,mean,set_printoptions
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  getms, writefermi, removezeros

from plotTools import plotxy
from pylab import *
from copy import deepcopy

def get_tetweights(dir):
    lines = nstrip(readfile(dir+'tetwgt'))
    print lines
    ntet = int(lines[0].split()[0])
    nbands = int(lines[0].split()[1])
#    nspin = int(lines[0][2]) #haven't added spin yet
    tetwgt = zeros((ntet,nbands),dtype = float64)
    for it in range(ntet):
        for ib in range(nbands):
            tetwgt[it,ib] = lines[it+1].split()[ib+1] #First line and first column are not weights             
    return ntet,nbands,tetwgt

def get_eigen(dir):
    lines = nstrip(readfile(dir+'eigenv'))
    nk = int(lines[0].split()[0])
    nbands = int(lines[0].split()[1])
#    nspin = int(lines[0][2]) E ignore for now  
    eigenv = zeros((nk,nbands),dtype = float64)
    for ik in range(nk):
        for ib in range(nbands):
            eigenv[ik,ib] = lines[ik+1].split()[ib+1]             
    return eigenv

def get_tetvecs(dir,ntet,nbands):
    lines = nstrip(readfile(dir+'tetvecs'))
    tetvecs = zeros((3,4,ntet),dtype = float64)
    for it in range(ntet):
        for ic in range(4):
            tetvecs[:,ic,it] = lines[5*it+ic+1].split()              
    return tetvecs

################# script #######################
dir = '/fslhome/bch/cluster_expansion/alal/test/f1_2/'

set_printoptions(precision=15)

os.chdir(dir)
ntet,nbands,tetwgt = get_tetweights(dir)
print ntet,tetwgt
ener = get_eigen(dir)
print ener
tvecs = get_tetvecs(dir,ntet,nbands)
print tvecs


    
#      
print 'Done'

