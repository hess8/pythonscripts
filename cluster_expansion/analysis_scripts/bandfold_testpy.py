#!/usr/bin/python
'''    
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  getms, writefermi, removezeros
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy
from pylab import *
from copy import deepcopy
fprec=float64

def read_eigenval(dir): 
    '''Read in k vectors and eigenvalues from vasp file'''
    os.chdir(dir)
    eigs = readfile('EIGENVAL')
    nb = int(eigs[5].split()[2])
    nk = int(eigs[5].split()[1])
    print nb, nk
    ks = zeros((nk,3))
    eners = zeros((nk,nb))   
    for ik in range(nk):
        istart = 7 + ik*(nb+2) 
        ks[ik,:]  = [float(eigs[istart].split()[0]), float(eigs[istart].split()[1]), float(eigs[istart].split()[2])]
        for ib in range(nb):
            eners[ik,ib] = float(eigs[istart+ib+1].split()[1])
    return [nb,nk,ks,eners]

title_detail =  'Cubic Al:Al (2.86 ang), cubic mesh, encut 500'
#title_detail =  'Cu:Cu, cubic mesh,f1-50,ediff 1e-7 '

dir1 = '/fslhome/bch/cluster_expansion/alal/cubic_al/equivk_c1-6_encut500/structs.cubmesh/c1_44/BANDS/'
dir2 = '/fslhome/bch/cluster_expansion/alal/cubic_al/equivk_c1-6_encut500/structs.cubmesh/c3_22/BANDS/'

[nb1,nk1,ks1,eners1] = read_eigenval(dir1)
[nb2,nk2,ks2,eners2] = read_eigenval(dir2)
eners1b = zeros((nk2,nb2))  #array with 2x the bands, 1/2 the kpoints 
#structure 1 is the one folded into the smaller BZ of the next structure.  
for ik in range(nk2):
    for ib in range(nb2/2):
        eners1b[ik,ib] = eners1[ik,ib] #copy normal ones to folded energy arrau
        eners1b[ik,ib + nb1] = eners1[ik + nk1/2,ib]
    eners1b[ik,:] = sorted(eners1b[ik,:])
print eners1b

#differences
diff = eners2 - eners1b
print diff

print eners2[0,:] ;print eners1b[0,:]

print 'Done'



