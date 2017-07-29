#!/usr/bin/python
'''
Comparison plot for different k mesh methods.  Allows reading external data from files for comparison
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean,sort,argsort,ceil,log10,int8,where
from numpy.linalg import norm
from analysisToolsVasp import getEnergy, getNkIBZ,getNatoms, readfile, writefile,electronicConvergeFinish
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts')
from symmetry import get_lattice_pointGroup, get_spaceGroup #these have vectors as ROWS
from kmeshroutines import readposcar

# from plotTools import plotxy
from pylab import figure,cm,plot,xlabel,ylabel,title,loglog,semilogy,legend,rcParams,style,rc
from matplotlib.colors import rgb2hex
# from matplotlib.colors import to_hex
from copy import deepcopy

def areEqual(x,y,eps):
    return abs(x-y)<eps

def plotData(datai,n,plotType,filter,doLegend,lablelStr):
    if plotType == 'linear':                      
        plot(datai['nKs'][:n],datai['eners'][:n],label=labelStr,\
              linestyle='None',color = datai['color'], marker = 'o',markeredgewidth=0.0)
    elif plotType == 'loglinear':                      
        semilogy(datai['nKs'][:n],datai['errs'][:n],label=labelStr,\
              linestyle='None',color = datai['color'], marker = 'o',markeredgewidth=0.0)
    elif plotType == 'loglog':                      
        loglog(datai['nKs'][:n],datai['errs'][:n],label=labelStr,\
              linestyle='None',color = datai['color'], marker = 'o',markeredgewidth=0.0)                             
    if doLegend:
        legend(loc='lower left',prop={'size':12});
        # show()
    fig.savefig('{}/{}_e_vs_n'.format(summaryPath,plotType))

testfile = 'POSCAR'

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/the99sym_newMethod',
#          '/fslhome/bch/cluster_expansion/vcmesh/vr_wc20']

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/vr_pow4',
#          '/fslhome/bch/cluster_expansion/vcmesh/vr_pow5',
#          '/fslhome/bch/cluster_expansion/vcmesh/vary_off05',
#          '/fslhome/bch/cluster_expansion/vcmesh/vr_pow7']

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/vary_off025',
#          '/fslhome/bch/cluster_expansion/vcmesh/vroff_04',
#          '/fslhome/bch/cluster_expansion/vcmesh/vroff_045',
#          '/fslhome/bch/cluster_expansion/vcmesh/vary_off05',
#          '/fslhome/bch/cluster_expansion/vcmesh/vary_off06',
#          '/fslhome/bch/cluster_expansion/vcmesh/vary_off075']

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/mt_fcc',
#          '/fslhome/bch/cluster_expansion/vcmesh/mt_bcc',
#           '/fslhome/bch/cluster_expansion/vcmesh/mt_cub',
#           '/fslhome/bch/cluster_expansion/vcmesh/mt_fo10']


paths = ['/fslhome/bch/cluster_expansion/vcmesh/mt_fcc']
# 
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/sc_fccOut',
#          '/fslhome/bch/cluster_expansion/vcmesh/sc_fcc']
# # r1
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/scond_fcc',
#          '/fslhome/bch/cluster_expansion/vcmesh/scond_fccOut',
#          '/fslhome/bch/cluster_expansion/vcmesh/sc_fccOut4',
#           '/fslhome/bch/cluster_expansion/vcmesh/sc_fccOutOF05',
#          '/fslhome/bch/cluster_expansion/vcmesh/sc_fOut03',
#          '/fslhome/bch/cluster_expansion/mpmesh/scond_mp']
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc04',
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc05',
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc06',
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc07',
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc08',
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc09', 
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc10', 
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc15',
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc20',   
#          '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc30',                           
#          '/fslhome/bch/cluster_expansion/mpmesh/scond_mp']

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/test','/fslhome/bch/cluster_expansion/mpmesh/semicond']

extpath = '/fslhome/bch/cluster_expansion/vcmesh/mueller_mp_data'
# extpath = None

coloring = 'method'
# coloring = 'indiv'
doLegend = True
doLabel = True

filter = '_' #string must be in dir name to be included
filter2 = None #'Cu_1' #for single structures.  set to None if using filter1 only
summaryPath = paths[0]

#count the number of plots:
iplot = 0
maxCalcs = 0
maxNk = 0
methods = []
for ipath,path in enumerate(paths):
    method = path.split('_')[-1]
    methods.append(method)
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
if not extpath is None:
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
                ('errs', '{}float'.format(maxCalcs)),('nKs', '{}int16'.format(maxCalcs)),('ns', '{}int8'.format(maxCalcs))])

# colorsList = []
# style.use('bmh')
# for i, item in enumerate(rcParams['axes.prop_cycle']):
#     colorsList.append(item['color']) 
style.use('fivethirtyeight')
# for i, item in enumerate(rcParams['axes.prop_cycle'][:-2]):
#     colorsList.append(item['color']) 

colorsList = [u'#30a2da', u'#fc4f30', u'#e5ae38', u'#6d904f', u'#8b8b8b',
              u'#348ABD', u'#A60628', u'#7A68A6', u'#467821', u'#D55E00', 
              u'#CC79A7', u'#56B4E9', u'#009E73', u'#F0E442', u'#0072B2']

colorsList = colorsList + ['b','m','y','c','k']
rcParams.update({'figure.autolayout': True})  
rcParams['axes.facecolor'] = 'white' 
rcParams['axes.linewidth'] = 1.0  
rcParams['axes.edgecolor'] = 'black' # axisbg=axescolor
rcParams['savefig.facecolor'] = 'white' # axisbg=axescolor
rcParams['lines.markersize'] = 4.5
#read all the data 
iplot = -1
for ipath, path in enumerate(paths): #my data
    print;print path
#     meshMethod = path.split('/')[-3][:3]+path.split('/')[-1][-3:]
    tag = path.split('/')[-1][-7:]
    os.chdir(path)
    if filter2 == None:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and filter in d])
    else:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d==filter2])
    nStructs = len(structs)
    print structs,path
    for istruct,struct in enumerate(structs):
#         print 'test', istruct, struct
        
        os.chdir(struct)
        if coloring == 'indiv':
#             if iplot < nplots -1:
            color = rgb2hex(cm.jet(1.*(iplot+1)/float(nplots)))
#             else:
#                 color = 'k' 
        elif coloring == 'method':
#             color =  colorsList[ipath]     
            color = None
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and os.path.exists('{}/OUTCAR'.format(d))])
        energies = []
        nKs = []
        nAtoms = []
        ns = [] #the base n of the run run
        nDone = 0
        for calc in calcs:  
            
            if electronicConvergeFinish(calc):
                ener = getEnergy(calc) #in energy/atom
                if not areEqual(ener,0,1e-5):
                    nDone +=1
                    energies.append(ener)
                    if 'vc' in path:
                        nK = getNkIBZ(calc,'KPOINTS')
                        
                    else:
                        nK = getNkIBZ(calc,'IBZKPT')
                    if nK > maxNk: maxNk = nK
                    nKs.append(nK)
                    ns.append(int(calc.split('_')[-1]))
        #sort by increasing number of kpoints
        if len(energies)>0: 
            iplot += 1
            nKs = array(nKs)
            energies = array(energies)
            ns = array(ns)
            order = argsort(nKs)
    #         print 'struct',struct
    #         print 'energies',energies
            energies = energies[order]
            ns = ns[order]
            nKs = sort(nKs)
            eref = energies[-1]#the last energy of each struct is that of the most kpoints
#             print 'struct',struct,'\t', mean(energies),'\t',eref,'\t', mean(energies-eref)*1000
    #         print 'eref',eref
    #         print 'energies sorted',energies, 'nKs', nKs
    #         print
           
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
            data[iplot]['ns'][:nDone] = ns
            data[iplot]['color'] = color
            method = path.split('_')[-1]
            data[iplot]['method'] = method
        os.chdir(path)
# os.chdir(extpath)
if not extpath is None:
    os.chdir(extpath)
    print; print atoms_methods
    for atom_method in atoms_methods:
        os.chdir(atom_method)
        if coloring == 'method':
            color = None
            if 'MP' in atom_method: 
#                 color = colorsList[len(paths)]
                method = 'MP'
                
            elif 'Mueller' in atom_method:
#                 color = colorsList[len(paths)+1]
                method = 'Mueller'
            if method not in methods:
                methods.append(method)
        if filter2 == None:
            structfiles = sorted([d for d in os.listdir(os.getcwd()) if os.path.getsize(d)>0])
        else:
            structfiles = sorted([d for d in os.listdir(os.getcwd()) if '_'.join(d.split('_')[:2])==filter2 and os.path.getsize(d)>0])
        for structfile in structfiles:
            if coloring == 'indiv':
                if iplot < nplots -1:
                    color = cm.jet(1.*(iplot+1)/float(nplots))
                else:
                    color = 'k'
            iplot += 1
            energies = []
            nKs = []
            lines = readfile(structfile)
            for line in lines:
                nK = int(line.split('\t')[0])
                if nK > maxNk: maxNk = nK
                nKs.append(nK)
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
nplots = iplot+1 

lines = [' ID , nKIBZ , ener , err, nAtoms, color\n']  
for iplot in range(nplots):
    n = data[iplot]['nDone']
    for icalc in range(n):#data[iplot]['eners'][:n].tolist()
        lines.append('{}_n{},{},{:15.12f},{:15.12f},{},{}\n'.format(data[iplot]['ID'], data[iplot]['ns'][icalc], data[iplot]['nKs'][icalc],\
             data[iplot]['eners'][icalc],data[iplot]['errs'][icalc],data[iplot]['nAtoms'],data[iplot]['color']))  
writefile(lines,'{}/summary.csv'.format(summaryPath)) 

#plots
if filter[0] == '_':filter = '' #labels can't begin with 
plotTypes = ['linear','loglinear','loglog']
ylabels = ['Vasp energy/atom (eV)','Error (meV)','Error (meV)']
xtext = 'N k-points'

for it,plotType in enumerate(plotTypes):
    fig = figure()
    ax1 = fig.add_subplot(111)
    xlabel(xtext)
    ylabel(ylabels[it]) 
    # title('Convergence vs mesh method')
    #ylim((1e-12,1e0))
    oldmethod = '' 
    for iplot in range(nplots):
        labelStr = None
        n = data[iplot]['nDone']
        if coloring == 'method':  
            method = data[iplot]['method'] 
            data[iplot]['color'] = colorsList[methods.index(method)] 
            if method != oldmethod:
                if doLabel: labelStr = '{} {}'.format(filter,data[iplot]['method'])
                plotData(data[iplot],n,plotType,filter,doLegend,labelStr)
                oldmethod = method;labelStr = None
            else:
                plotData(data[iplot],n,plotType,filter,doLegend,labelStr)
        elif coloring == 'indiv': 
            if doLabel: labelStr = '{} {}'.format(filter,data[iplot]['ID'])
            plotData(data[iplot],n,plotType,filter,doLegend,labelStr)
 
#Method averaging
if coloring == 'method':
    smoothFactor = 4.0
    print 'Averaging, plotting method errors'
    fig = figure()
    ax1 = fig.add_subplot(111)
    xlabel('N k-points (smoothed by factor {})'.format(int(smoothFactor)))
    ylabel('Error (meV)') 
    for im,method in enumerate(methods):
        methnKmax = 0 
        nbins = int(10*ceil(log10(maxNk)))# 10 bins per decade
        nKbins = array([(10.0**(1/10.0))**i for i in range(nbins)])
        binCounts = zeros(nbins,dtype = int8)
        binErrs = zeros(nbins,dtype = float)       
        for iplot in range(nplots):
            if data[iplot]['method'] == method:
                for icalc in range(data[iplot]['nDone']-1):
                    nK = data[iplot]['nKs'][icalc]
                    if nK>methnKmax: methnKmax = nK
                    if nK>1:
                        for ibin in range(nbins):
                            if abs(log10(nK/nKbins[ibin])) <= log10(smoothFactor)\
                               and nKbins[ibin]<= maxNk:
                               binErrs[ibin] += data[iplot]['errs'][icalc]
                               binCounts[ibin] += 1
        mask = where(binCounts>0)
        binErrs = binErrs[mask[0]]
        binCounts = binCounts[mask[0]]
        nKbins = nKbins[mask[0]]
        nbins = len(nKbins)
        avgErrs = [binErrs[ibin]/binCounts[ibin] for ibin in range(nbins)]
        #plot
        loglog(nKbins,avgErrs,label = method,\
              color = colorsList[im], marker = None)
        print 'Method',method, 'nKmax',methnKmax
    legend(loc='lower left',prop={'size':12});
    fig.savefig('{}/methodErrs'.format(summaryPath))
        
        #print 'Method {} \terror {}'.format(method,normErr), 'counts',mcounts[im]
print 'Done'

