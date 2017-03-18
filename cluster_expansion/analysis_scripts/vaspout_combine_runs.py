#!/usr/bin/python
'''
Comparison plot for different k mesh methods
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  writefermi, removezeros
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy
from pylab import *
from copy import deepcopy
fprec=float64

testfile = 'POSCAR'

def getibest(dirs):
#    mrange = []
    mmax = 0
    for i,dir in enumerate(dirs):
        if dir[1]=='1' and dir[2]=='_':
            m = int(dir.split('_')[-1])
#            print dir, m
#            if m > mmax:
#                ibest = i
#                mmax = m
            if m == 24:
                ibest = i
                mmax = m
#            mrange.append(m)
    return mmax, ibest

################# script #######################
################# script #######################v
################# script #######################

paths = ['/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/AFLOWDATAn/Cu_pvPt']
summaryPath = ['/fslhome/bch/cluster_expansion/vcmesh/']
iplot = 0
maxCalcs = 0
#count the number of total runs:
for ipath,maindir in enumerate(paths):
    os.chdir(maindir)
    structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
    for struct in struct:
        iplot += 1
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        if len(calcs)>maxCalcs: maxCalcs = len(calcs)
        for calc in calcs:
            ncalcs += 1   
nplots = iplot       
data = zeros((nplots,),dtype = [('ID', 'S15'),('nDone','int8'),('eners', '{}*float'.format(maxCalcs)),\
                ('errs', '{}*float'.format(maxCalcs)),('nKs', '{}*int8'.format(maxCalcs))])
#read all the data            
iplot = 0
for ipath ,maindir in enumerate(paths):
    titles.append(maindir.split('/')[-3][0].upper()+path.split('/')[-3][1].lower())
    meshMethod = maindirsplit('/')[-4]
    os.chdir(maindir)
    structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
    nStructs = len(structs)
    print structs
    for struct in struct:
        iplot += 1
        os.chdir(struct)
        runs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        energies = []
        nKs = []
        nDone = 0
        for run in runs:  
            if finished(run):
                nDone +=1
                ener = getEnergy(run) #in energy/atom
                if notEquals(ener,0):
                    energies.append(ener)
                    nKs.append(getNk(run))
        energies = array(energies)
        eref = energies[-1] #the last energy of each struct should be the most kpoints
        errs = (energies-eref)*1000 + 1e-16 #now in meV, with 
        plotErrs()
        data[iplot]['ID'] = '{} {}'.format(struct,meshMethod)
        data[iplot]['nDone'] = nDone
        data[iplot]['eners'][:len(runs)] = energies
        data[iplot]['errs'][:len(runs)] = errs
        data[iplot]['nKs'][:len(runs)] = nKs           
        os.chdir(maindir)
    
#     #log(err) vs NkIBZ
#     fig = figure()
#     semilogy(NkIBZ,err,'ro')
#     title(titleadd + ' Error vs Nk in IBZKPT')
#     xlabel('Nk')
#     ylabel('error')   
#     fig.savefig('nk_log_err')  
#     
#     #log(err) vs log(NkIBZ)
#     fig = figure()
#     loglog(NkIBZ,err,'ro')
#     title(titleadd + ' Error vs Nk in IBZKPT')
#     xlabel('Nk')
#     ylabel('error')   
#     fig.savefig('nk_loglog_err') 

fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('N kpoints')
ylabel('Error (meV)') 
title('Convergence vs mesh method')

#ylim((1e-12,1e0))
for iplot in range(nplots):
    n = data[iplot]['nDone']  
    ax1.semilogy(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
      label=data[iplot]['ID'],linestyle='None',color=cm.jet(1.*(iplot+1)/float(nplots)), marker = 'o') 
plt.legend(loc='upper right',prop={'size':14});
show()
fig.savefig('{}/log_err_vs_n'.format(summaryPath,string)) 
print 'Done'

