#!/usr/bin/python
'''
Comparison plot for different k mesh methods
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean,sort,argsort
from numpy.linalg import norm
from analysisToolsVasp import getEnergy, getNkIBZ,getNatoms, readfile, writefile,electronicConvergeFinish
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts')
from symmetry import get_lattice_pointGroup, get_spaceGroup #these have vectors as ROWS
from kmeshroutines import readposcar

# from plotTools import plotxy
from pylab import *
from copy import deepcopy

def areEqual(x,y,eps):
    return abs(x-y)<eps

testfile = 'POSCAR'


# 
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/limNvary_off05',
#          '/fslhome/bch/cluster_expansion/vcmesh/limNvary_off-05']
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/the99sym','/fslhome/bch/cluster_expansion/mpmesh/mpPure']
paths = ['/fslhome/bch/cluster_expansion/vcmesh/scond','/fslhome/bch/cluster_expansion/mpmesh/semicond']
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/test','/fslhome/bch/cluster_expansion/mpmesh/semicond']

filter = 'C' #string must be in dir name to be included
summaryPath = paths[0]
# summaryPath = '/fslhome/bch/cluster_expansion/vcmesh/cu17Jul17/'
# summaryPath = paths[1]
iplot = 0
maxCalcs = 0
#count the number of plots:
for ipath,maindir in enumerate(paths):
    os.chdir(maindir)
    structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and filter in d])
#     structs = ['f9292']; print'only struct 9292'
    for struct in structs:
        os.chdir(struct)
        iplot += 1
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        if len(calcs)>maxCalcs: maxCalcs = len(calcs)
        os.chdir(maindir)
nplots = iplot 
if nplots < len(paths): sys.exit('Stop.  Not enough structures match filter')      
data = zeros(nplots,dtype = [('ID', 'S15'),('nDone','int8'),('nAtoms','int8'),('nops','int8'),('eners', '{}float'.format(maxCalcs)),\
                ('errs', '{}float'.format(maxCalcs)),('nKs', '{}int16'.format(maxCalcs))])
#read all the data 
iplot = -1
for ipath, maindir in enumerate(paths):
#     meshMethod = maindir.split('/')[-3][:3]+maindir.split('/')[-1][-3:]
    meshMethod = maindir.split('/')[-1][-7:]
    os.chdir(maindir)
    structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and filter in d])
#     structs = ['f9292']; print'only struct 9292'
    nStructs = len(structs)
    print;print structs,maindir
    for istruct,struct in enumerate(structs):
#         print 'test', istruct, struct
        iplot += 1
        os.chdir(struct)
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and os.path.exists('{}/OUTCAR'.format(d))])
        energies = []
        nKs = []
        nAtoms = []
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
        if len(energies)>0: 
            nKs = array(nKs)
            energies = array(energies)
            order = argsort(nKs)
    #         print 'struct',struct
    #         print 'energies',energies
            energies = energies[order]
            eref = energies[-1]#the last energy of each struct is that of the most kpoints
    #         print 'eref',eref
    #         print 'energies sorted',energies, 'nKs', nKs
    #         print
            nKs = sort(nKs)
    #         print 'NKs', nKs      
    #         if ipath == 0 and istruct == 0: ******* useful if comparing methods on the same struct with a trusted first method*******
    #             eref = energies[-1] #the last energy of each struct is that of the most kpoints
    
            errs = abs(energies-eref)*1000 + 1e-6 #now in meV 
    #         plotErrs()
            data[iplot]['ID'] = '{} {}'.format(struct,meshMethod)
            nAtoms = getNatoms('{}/POSCAR'.format(calc))
            data[iplot]['nAtoms'] = nAtoms
            data[iplot]['nDone'] = nDone
            data[iplot]['eners'][:nDone] = energies
            data[iplot]['errs'][:nDone] = errs
            data[iplot]['nKs'][:nDone] = nKs
        os.chdir(maindir)
lines = [' ID , nKIBZ , ener , err, nAtoms \n']  
for iplot in range(nplots):
    n = data[iplot]['nDone']
    for icalc in range(n):#data[iplot]['eners'][:n].tolist()
        lines.append('{},{},{:15.12f},{:15.12f},{}\n'.format(data[iplot]['ID'], data[iplot]['nKs'][icalc],\
             data[iplot]['eners'][icalc],data[iplot]['errs'][icalc],data[iplot]['nAtoms']))  
writefile(lines,'{}/summary.csv'.format(summaryPath))   
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('N kpoints')
ylabel('Vasp energy/atom (eV)') 
title('Convergence vs mesh method')
#ylim((1e-12,1e0))
for iplot in range(nplots):
    if iplot < nplots -1:
        plotcolor = cm.jet(1.*(iplot+1)/float(nplots))
    else:
        plotcolor = 'k'
    n = data[iplot]['nDone']  
#     print 'iplot',data[iplot]['eners'][:n], data[iplot]['nKs'][:n]
    plot(data[iplot]['nKs'][:n],data[iplot]['eners'][:n],\
      label=data[iplot]['ID'],linestyle='None',color = plotcolor, marker = 'o',markeredgewidth=0.0) 
legend(loc='lower left',prop={'size':12});
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
    if iplot < nplots -1:
        plotcolor = cm.jet(1.*(iplot+1)/float(nplots))
    else:
        plotcolor = 'k'
    n = data[iplot]['nDone']  
    semilogy(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
      label=data[iplot]['ID'],linestyle='None',color=plotcolor, marker = 'o',markeredgewidth=0.0) 
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
    if iplot < nplots -1:
        plotcolor = cm.jet(1.*(iplot+1)/float(nplots))
    else:
        plotcolor = 'k'
    n = data[iplot]['nDone']  
    ax1.loglog(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
      label=data[iplot]['ID'],linestyle='None',color=plotcolor, marker = 'o',markeredgewidth=0.0) 
legend(loc='lower left',prop={'size':12});
# show()
fig.savefig('{}/loglog_err_vs_n'.format(summaryPath)) 


print 'Done'

