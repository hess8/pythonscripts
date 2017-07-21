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
from pylab import figure,cm,plot,xlabel,ylabel,title,loglog,semilogy,legend,rcParams,style,rc
from copy import deepcopy

def areEqual(x,y,eps):
    return abs(x-y)<eps

testfile = 'POSCAR'

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/the99sym_newMethod']
paths = ['/fslhome/bch/cluster_expansion/vcmesh/vary_off05',
         '/fslhome/bch/cluster_expansion/vcmesh/vary_off00',
         '/fslhome/bch/cluster_expansion/vcmesh/vary_off-05']

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/the99sym_newMethod','/fslhome/bch/cluster_expansion/mpmesh/mpPure_MPbch']
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/semicond','/fslhome/bch/cluster_expansion/mpmesh/semicond']
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/test','/fslhome/bch/cluster_expansion/mpmesh/semicond']
extpath = '/fslhome/bch/cluster_expansion/vcmesh/mueller_mp_data'

filter = '_' #string must be in dir name to be included
filter2 = None #'Cu_1' #for single structures.  set to None if using filter1 only
summaryPath = paths[0]
# summaryPath = '/fslhome/bch/cluster_expansion/vcmesh/cu17Jul17/'
# summaryPath = paths[1]

#count the number of plots:
iplot = 0
maxCalcs = 0
for ipath,path in enumerate(paths):
    os.chdir(path)
    if filter2 == None:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and filter in d])
    else:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d==filter2])
    for struct in structs:
        os.chdir(struct)
        iplot += 1
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        if len(calcs)>maxCalcs: maxCalcs = len(calcs)
        os.chdir(path)

#external data is of the form extpath/atom_method/struct.csv.  The csv has energies vs nK
os.chdir(extpath)
atoms_methods = sorted([d for d in os.listdir(extpath) if os.path.isdir(d) and filter in d])# os.chdir(extpath)
for atom_method in atoms_methods:
    atom = atom_method.split('_')[0]
    os.chdir(atom_method)
    os.system('rm -r .*lock*')
    for structfile in os.listdir(os.getcwd()):
        if atom not in structfile:
            os.system('mv {} {}_{}'.format(structfile,atom,structfile)) #so that file has atom name at beginning
    if filter2 == None:
        structfiles = sorted([d for d in os.listdir(os.getcwd()) if os.path.getsize(d)>0])
    else:
        structfiles = sorted([d for d in os.listdir(os.getcwd()) if '_'.join(d.split('_')[:2])==filter2 and os.path.getsize(d)>0])
    for structfile in structfiles:
        iplot += 1
        #count number of points in this structfile
        lines = readfile(structfile)
        if len(lines)>maxCalcs: maxCalcs = len(lines)     
    os.chdir(extpath)

nplots = iplot 
if nplots < len(paths): sys.exit('Stop.  Not enough structures match filter')      
data = zeros(nplots,dtype = [('ID', 'S15'),('color', 'S15'),('method', 'S15'),('nDone','int8'),('nAtoms','int8'),('nops','int8'),('eners', '{}float'.format(maxCalcs)),\
                ('errs', '{}float'.format(maxCalcs)),('nKs', '{}int16'.format(maxCalcs))])
style.use('fivethirtyeight')
colors = []
for i, item in enumerate(rcParams['axes.prop_cycle']):
    colors.append(item['color']) 
rcParams.update({'figure.autolayout': True})  
rcParams['axes.facecolor'] = 'white' 
rcParams['axes.linewidth'] = 1.0  
rcParams['axes.edgecolor'] = 'black' # axisbg=axescolor
rcParams['savefig.facecolor'] = 'white' # axisbg=axescolor
#read all the data 
iplot = -1
os.chdir(extpath)
print; print atoms_methods
for atom_method in atoms_methods:
    os.chdir(atom_method)
    if 'MP' in atom_method: 
        color = colors[len(paths)]
        method = 'MP'
    elif 'Mueller' in atom_method:
        color = colors[len(paths)+1]
        method = 'Mueller'
    if filter2 == None:
        structfiles = sorted([d for d in os.listdir(os.getcwd()) if os.path.getsize(d)>0])
    else:
        structfiles = sorted([d for d in os.listdir(os.getcwd()) if '_'.join(d.split('_')[:2])==filter2 and os.path.getsize(d)>0])
    for structfile in structfiles:
#         print 'file',structfile
        iplot += 1
        energies = []
        nKs = []
        lines = readfile(structfile)
        for line in lines:
            nKs.append(int(line.split('\t')[0]))
            energies.append(-float(line.split('\t')[1].split('\r')[0]))
        nKs = array(nKs)
        energies = array(energies)
        nDone = len(energies)
        order = argsort(nKs)
        energies = energies[order]
        eref = energies[-1]#the last energy of each struct is that of the most kpoints
        nKs = sort(nKs)
        errs = abs(energies-eref)*1000 + 1e-6 #now in meV 
        struct = structfile.split('_')[0]
        data[iplot]['ID'] = atom_method + struct
        data[iplot]['nAtoms'] = 0
        data[iplot]['nDone'] = len(energies)
        data[iplot]['eners'][:nDone] = energies
        data[iplot]['errs'][:nDone] = errs
        data[iplot]['nKs'][:nDone] = nKs
        data[iplot]['color'] = color
        data[iplot]['method'] = method
    os.chdir(extpath)
for ipath, path in enumerate(paths):
#     meshMethod = path.split('/')[-3][:3]+path.split('/')[-1][-3:]
    tag = path.split('/')[-1][-7:]
    os.chdir(path)
    if filter2 == None:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and filter in d])
    else:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d==filter2])
    nStructs = len(structs)
    print;print structs,path
    for istruct,struct in enumerate(structs):
#         print 'test', istruct, struct
        
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
                    if 'vc' in path:
                        nK = getNkIBZ(calc,'KPOINTS')
                    else:
                        nK = getNkIBZ(calc,'IBZKPT')
                    nKs.append(nK)
        #sort by increasing number of kpoints
        if len(energies)>0: 
            iplot += 1
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
            data[iplot]['ID'] = '{} {}'.format(struct,tag)
            nAtoms = getNatoms('{}/POSCAR'.format(calc))
            data[iplot]['nAtoms'] = nAtoms
            data[iplot]['nDone'] = nDone
            data[iplot]['eners'][:nDone] = energies
            data[iplot]['errs'][:nDone] = errs
            data[iplot]['nKs'][:nDone] = nKs
            data[iplot]['color'] = colors[ipath]
            method = path.split('_')[-1]
            data[iplot]['method'] = method
        os.chdir(path)
# os.chdir(extpath)

nplots = iplot+1 

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
# title('Convergence vs mesh method')
#ylim((1e-12,1e0))
oldmethod = ''
methods = []
if filter[0] == '_':filter = '' #labels can't begin with _
for iplot in range(nplots):
    n = data[iplot]['nDone']     
    method = data[iplot]['method']
    if method != oldmethod and method not in methods:
        methods.append(method)
        plot(data[iplot]['nKs'][:n],data[iplot]['eners'][:n],\
          label='{} {}'.format(filter,data[iplot]['method']),linestyle='None',color = data[iplot]['color'], marker = 'o',markeredgewidth=0.0)
    else:
        plot(data[iplot]['nKs'][:n],data[iplot]['eners'][:n],\
          linestyle='None',color = data[iplot]['color'], marker = 'o',markeredgewidth=0.0)     
legend(loc='lower left',prop={'size':12});
# show()
fig.savefig('{}/energy_vs_n'.format(summaryPath))         

oldmethod = ''
methods = []        
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('N kpoints')
ylabel('Error (meV)') 
# title('Convergence vs mesh method')
#ylim((1e-12,1e0))
for iplot in range(nplots):
    n = data[iplot]['nDone'] 
    method = data[iplot]['method']
    if method != oldmethod and method not in methods:
        methods.append(method)
        semilogy(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
          label='{} {}'.format(filter,data[iplot]['method']),linestyle='None',color = data[iplot]['color'], marker = 'o',markeredgewidth=0.0)
    else:
        semilogy(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
          linestyle='None',color = data[iplot]['color'], marker = 'o',markeredgewidth=0.0)     
legend(loc='lower left',prop={'size':12});
# show()
fig.savefig('{}/log_err_vs_n'.format(summaryPath)) 

#log-log
oldmethod = ''
methods = []
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('N kpoints')
ylabel('Error (meV)') 
# title('Convergence vs mesh method')
for iplot in range(nplots):
    n = data[iplot]['nDone']
    method = data[iplot]['method']
    if method != oldmethod and method not in methods:
        methods.append(method)
        loglog(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
          label='{} {}'.format(filter,data[iplot]['method']),linestyle='None',color = data[iplot]['color'], marker = 'o',markeredgewidth=0.0)
    else:
        loglog(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
          linestyle='None',color = data[iplot]['color'], marker = 'o',markeredgewidth=0.0)   


  
#     ax1.loglog(data[iplot]['nKs'][:n],data[iplot]['errs'][:n],\
#       label=data[iplot]['ID'],linestyle='None',color=plotcolor, marker = 'o',markeredgewidth=0.0) 
    
    
legend(loc='lower left',prop={'size':12});
# show()
fig.savefig('{}/loglog_err_vs_n'.format(summaryPath)) 


print 'Done'

