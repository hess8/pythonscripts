#!/usr/bin/python
'''
Comparison plot for different k mesh methods
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean,sort,argsort
from numpy.linalg import norm
from analysisToolsVasp import getEnergy, getNkIBZ, readfile, writefile,electronicConvergeFinish
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
# from plotTools import plotxy
from pylab import *
from copy import deepcopy

def areEqual(x,y,eps):
    return abs(x-y)<eps

testfile = 'POSCAR'

# paths = ['/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/AFLOWDATAn/Cu_pvPt',\
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/AFLOWDATAn/Cu_pvPt']
# paths = ['/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/cubicTest',\
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTest']
#  
# paths = ['/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/AFLOWDATAn/Cu_pvPt',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestBCC',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestCUB',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestFCC'
#          ]
# 
# paths = ['/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/cubicTest',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestNoRedistrMv',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestNoRedistr',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrFCC',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrBCC',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrCUB',
#          '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestNoSymComm']

paths = ['/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/f1MP',
         '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1DP']


summaryPath = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/'
# summaryPath = paths[1]
iplot = 0
maxCalcs = 0
#count the number of plots:
for ipath,maindir in enumerate(paths):
    os.chdir(maindir)
    structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
#     structs = ['f9292']; print'only struct 9292'
    for struct in structs:
        os.chdir(struct)
        iplot += 1
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        if len(calcs)>maxCalcs: maxCalcs = len(calcs)
        os.chdir(maindir)
nplots = iplot       
data = zeros(nplots,dtype = [('ID', 'S15'),('nDone','int8'),('eners', '{}float'.format(maxCalcs)),\
                ('errs', '{}float'.format(maxCalcs)),('nKs', '{}int16'.format(maxCalcs))])
#read all the data            
iplot = -1
for ipath, maindir in enumerate(paths):
#     meshMethod = maindir.split('/')[-3][:3]+maindir.split('/')[-1][-3:]
    meshMethod = maindir.split('/')[-1][-7:]
    os.chdir(maindir)
    structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
#     structs = ['f9292']; print'only struct 9292'
    nStructs = len(structs)
    print;print structs,maindir
    for istruct,struct in enumerate(structs):
        iplot += 1
        os.chdir(struct)
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        energies = []
        nKs = []
        nDone = 0
        for calc in calcs:  
#             print 'calc',calc
            if electronicConvergeFinish(calc):
                ener = getEnergy(calc) #in energy/atom
                if not areEqual(ener,0,1e-5):
                    nDone +=1
                    energies.append(ener)
                    if 'vc' in maindir:
                        nK = getNkIBZ(calc,'KPOINTS')
                    else:
                        nK = getNkIBZ(calc,'IBZKPT')
                    nKs.append(nK)
        #sort by increasing number of kpoints
        nKs = array(nKs)
        energies = array(energies)
        order = argsort(nKs)
        energies = energies[order]
        nKs = sort(nKs)
#         print 'energies',  energies
#         print 'NKs', nKs      
#         eref = energies[-1] #the last energy of each struct is that of the most kpoints
        if ipath == 0 and istruct == 0:
            eref = energies[-1] #the last energy of each struct is that of the most kpoints
#         eref = energies[-1] #the last energy of each struct is that of the most kpoints
        errs = abs(energies-eref)*1000 + 1e-6 #now in meV 
#         plotErrs()
        data[iplot]['ID'] = '{} {}'.format(struct,meshMethod)
        data[iplot]['nDone'] = nDone
        data[iplot]['eners'][:nDone] = energies
        data[iplot]['errs'][:nDone] = errs
        data[iplot]['nKs'][:nDone] = nKs
        os.chdir(maindir)
lines = [' ID , nKIBZ , ener , err \n']  
for iplot in range(nplots):
    n = data[iplot]['nDone']
    for icalc in range(n):#data[iplot]['eners'][:n].tolist()
        lines.append('{},{},{:15.12f},{:15.12f}\n'.format(data[iplot]['ID'], data[iplot]['nKs'][icalc],\
             data[iplot]['eners'][icalc],data[iplot]['errs'][icalc]))   
writefile(lines,'{}/summary.csv'.format(summaryPath))   
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('N kpoints')
ylabel('Vasp energy/atom (eV)') 
title('Convergence vs mesh method')
#ylim((1e-12,1e0))
for iplot in range(nplots):
    n = data[iplot]['nDone']  
#     print 'iplot',data[iplot]['eners'][:n], data[iplot]['nKs'][:n]
    plot(data[iplot]['nKs'][:n],data[iplot]['eners'][:n],\
      label=data[iplot]['ID'],linestyle='None',color=cm.jet(1.*(iplot+1)/float(nplots+1)), marker = 'o',markeredgewidth=0.0) 
legend(loc='upper center',prop={'size':12});
# show()
fig.savefig('{}/energy_vs_n'.format(summaryPath))         
        
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('N kpoints')
ylabel('Error (meV)') 
title('Convergence vs mesh method')
#ylim((1e-12,1e0))
for iplot in range(nplots):
    n = data[iplot]['nDone']  
    semilogy(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
      label=data[iplot]['ID'],linestyle='None',color=cm.jet(1.*(iplot+1)/float(nplots+1)), marker = 'o',markeredgewidth=0.0) 
legend(loc='lower left',prop={'size':12});
# show()
fig.savefig('{}/log_err_vs_n'.format(summaryPath)) 

#log-log
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('N kpoints')
ylabel('Error (meV)') 
title('Convergence vs mesh method')
#ylim((1e-12,1e0))
for iplot in range(nplots):
    n = data[iplot]['nDone']  
    ax1.loglog(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
      label=data[iplot]['ID'],linestyle='None',color=cm.jet(1.*(iplot+1)/float(nplots+1)), marker = 'o',markeredgewidth=0.0) 
legend(loc='lower left',prop={'size':12});
# show()
fig.savefig('{}/loglog_err_vs_n'.format(summaryPath)) 


print 'Done'

