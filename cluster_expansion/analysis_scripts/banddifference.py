#!/usr/bin/python
'''   
Compares two sets of band structure (identical number of bands and kpoints), and denotes the differences by color of the bands
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  getms, writefermi, removezeros,getEf
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy, colorline
from matplotlib.pyplot import pcolor, colorbar
from pylab import *
from copy import deepcopy
fprec=float64

def getEf(folder): 
    '''Finds fermi energy from OUTCAR'''
    lastdir = os.getcwd()
    os.chdir(folder)
    try:
        outcar = open('OUTCAR','r')
        text = outcar.readlines()
        proc = subprocess.Popen(['grep','-i','E-fermi','OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        ef = newstring[0].split()[2]
        print folder, ef
    except:
        ef = str(0.00)
    os.chdir(lastdir)
    return ef

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

title_detail =  'Cubic Al:Al (2.86 ang), c3 (4x4x2) when stretched 0.5% along y and z) \n Lowest 5 bands. Color: log10(diff(eV)+1e-6)'
plotfile = 'bands_error_c3_c3stretched'
#title_detail =  'Cu:Cu, cubic mesh,f1-50,ediff 1e-6 '

dir1 = '/fslhome/bch/cluster_expansion/alal/cubic_al/equivk_c1-6_encut500/structs.cubmesh/c3_22/BANDS/'
dir2 = '/fslhome/bch/cluster_expansion/alal/cubic_al/equivk_c1-6_encut500/structs.cubmesh/c3_22/stretched/BANDS/'

ef = getEf(dir1)

[nb1,nk1,ks1,eners1] = read_eigenval(dir1)
[nb2,nk2,ks2,eners2] = read_eigenval(dir2)
if nb1 != nb2:
    print nb1, nb2
    sys.exit('Number of bands is not the same')
elif nk1 != nk2:
    print nk1, nk2
    sys.exit('Number of kpoints is not the same')    

#differences
diff = abs(eners1 - eners2)
print diff

print eners1[0,:] ;print eners2[0,:]
diffmin = amin(diff[:,0:nb1/2]), amax(diff[:,0:nb1/2]) #lower half of bands

#plot only the lower half of the bands...the highes energy ones are not occupied, and have strange differences

#flatten the data so we can do a scatter plot
X = zeros((nk1*nb1/2,1))
Y = zeros((nk1*nb1/2,1))
Z = zeros((nk1*nb1/2,1))
for ik in range(nk1):
    for ib in range(nb1/2):
        i = ik*nb1/2 + ib #where we are in the bands line list
        X[i] = ik
        Y[i] = eners2[ik,ib]
        Z[i] = log10(diff[ik,ib] + 1e-6)
        
print 'X',X[0:50]
print 'Y',Y[0:50]
print 'Z',Z[0:50]
print amin(Z), amax(Z)

ticklabels = [r'\Gamma', 'X', 'M', r'\Gamma', 'Z','R','A','Z X','R M','A']
nsegments = (len(ticklabels)-1)
k_per_segment =  nk1/nsegments
ticks = [i*k_per_segment for i in range(nsegments+1)]

print ticks

fig, axes = subplots()
scatter(X,Y,c=Z, cmap='gnuplot',vmin=amin(Z), vmax=amax(Z), edgecolor='None',)
plot([0,nk1],[ef,ef])
xlim((0,nk1))
ylim((amin(Y)-1,amax(Y)+1))
for x in ticks: #lines for divisions
    plot([x,x],[amin(Y)-1,amax(Y)+1],color = 'black')
xticks(ticks, ['$%s$' % n for n in [r'\Gamma', 'X', 'M', r'\Gamma', 'Z','R','A','Z X','R M','A']])
ylabel('energy (eV)')
colorbar()
title(title_detail)
show() 
os.chdir(dir1)
fig.savefig(plotfile)

print 'Done'



