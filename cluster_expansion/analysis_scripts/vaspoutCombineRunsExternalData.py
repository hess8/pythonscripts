#!/usr/bin/python
'''
Comparison plot for different k mesh methods.  Allows reading external data from files for comparison
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean,sort,argsort,ceil,log10,int8,int32,where
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

# def writeSym(self):
#         writefile(['nops: {}\n'.format(self.nops),'IBZvolCut: {}\n'.format(self.IBZvolCut)],'sym.out')
def readSym(dir):
    lines = readfile('{}/sym.out'.format(dir))
    nops = int(lines[0].split(':')[1])
    IBZvolCut = float(lines[1].split(':')[1])
#     IBZvol = float(lines[2].split(':')[1])
    IBZvol = None
    return nops, IBZvolCut, IBZvol

def copyData(struct,data):
    n = struct.split('_')[1]
    for i,ID in enumerate(data['ID']):
        if struct in ID and n in ID:
            icopy = i
            break
    return data[icopy]['nops'],data[icopy]['IBZvolcut'],data[icopy]['nAtoms']

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
#         loglog(datai['nKs'][:n],datai['errs'][:n],label=labelStr,\
#               linestyle='-',color = datai['color'], marker = 'o',markeredgewidth=0.0)                              
    if doLegend:
        legend(loc='lower left',prop={'size':12});
        # show()
    fig.savefig('{}/{}_e_vs_n'.format(summaryPath,plotType))

testfile = 'POSCAR'

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_mxwc10',
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_mxwc25',
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_testfcc',
#          '/fslhome/bch/cluster_expansion/mpmesh/scond_mp']

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

paths = [ 
# '/fslhome/bch/cluster_expansion/vcmesh/mt_fcc',
        '/fslhome/bch/cluster_expansion/vcmesh/mt_grid17Aug/bestRun' ]
# ,
#           '/fslhome/bch/cluster_expansion/vcmesh/mt_cub',
#           '/fslhome/bch/cluster_expansion/vcmesh/mt_fo10']

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_testfcc',
#     '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_fcc',
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_bcc',
#           '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_cub']

# paths = ['/fslhome/bch/cluster_expansion/vcmesh/mt_grid/r0']
    
#     '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_testfccParams', 
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_bcc', 
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_cub', 
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_junk',
#          '/fslhome/bch/cluster_expansion/mpmesh/scond_mp']
#
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_fccOut',
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_fcc']
# # r1
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/semiconductors/scond_fcc',
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/scond_fccOut',
#          '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_fccOut4',
#           '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_fccOutOF05',
#          '/fslhome/bch/cluster_expansion/vcmesh/sc_fOut03',
#          '/fslhome/bch/cluster_expansion/mpmesh/scond_mp']

# paths = ['/fslhome/bch/cluster_expransion/vcmesh/test','/fslhome/bch/cluster_expansion/mpmesh/semicond']

extpaths = ['/bluehome/bch/fsl_groups/fslg_datamining/Mueller']
# extpath = None
useSym = False
coloring = 'method'
# coloring = 'indiv'
doLegend = True
doLabel = True
collateMeshMat = True  #gather the mathematica plotting code for the spheres and IBZ for each calculation in a file for each method. 
smoothFactor = 2.0
filter = 'Al_' #string must be in dir name to be included
filter2 = None #'Cu_1' #for single structures.  set to None if using filter1 only
summaryPath = paths[0]

#count the number of plots:
iplot = 0
maxCalcs = 0
maxNk = 0
methods = []
for ipath,path in enumerate(paths):
    method = path.split('_')[-1].split('/')[0]
    methods.append(method)
    os.chdir(path)
    if collateMeshMat:
        meshPlots = open('IBZmeshPlots','w')
    if filter2 == None:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and filter in d])
    else:
        structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d==filter2])
    for struct in structs:
        if collateMeshMat:meshPlots.write('(* {} *)\n'.format(struct))
        os.chdir(struct)
        iplot += 1
        calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        if len(calcs) > maxCalcs: maxCalcs = len(calcs)
        if collateMeshMat:
            for calc in calcs:
                os.chdir(calc)
                if os.path.exists('cell_IBZmesh.m'):
                    meshPlots.write('\t(* {} *)\n'.format(calc))
                    lines = readfile('cell_IBZmesh.m')
                    lines.append( '\n\n')
                    meshPlots.writelines(lines)
                os.chdir('../')          
        os.chdir(path)
    if collateMeshMat:
        meshPlots.close()
        
#external run paths are of the form extmethodpath/atom_convergence/11_atom 
if not extpaths is None:
    for ipath,extpath in enumerate(extpaths):
        os.chdir(extpath)
        method = extpath.split('/')[0]
        atomdirs = sorted([d for d in os.listdir(extpath) if os.path.isdir(d) and filter in d])# os.chdir(extpath)
        for dir in atomdirs:
            atom = dir.split('_')[0]
            os.chdir(dir)
            structs = []
            for item in os.listdir(os.getcwd()): 
                if os.path.isdir(item):
                    if item[0].isdigit():
                        os.system('mv {} {}_{}'.format(item,atom,item.split('_')[0])) #rename to format 'Cu_11' instead of 11_atom
            for item in os.listdir(os.getcwd()): 
                if os.path.isdir(item) and (filter2 == None or filter2 in struct):
                    structs.append(item)
            for struct in structs:
                os.chdir(struct)
                iplot += 1
                calcs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
                if len(calcs)>maxCalcs: maxCalcs = len(calcs)        
                os.chdir('../')
            os.chdir(extpath)               

nplots = iplot 
if nplots < len(paths): sys.exit('Stop.  Structures do not match filter')      
data = zeros(nplots,dtype = [('ID', 'S25'),('color', 'S15'),('method', 'S15'),\
                             ('nDone','int32'),('nAtoms','int32'),('nops','int8'),\
                             ('IBZvolcut','float'),('IBZvol','float'),\
                             ('eners', '{}float'.format(maxCalcs)), ('errs', '{}float'.format(maxCalcs)),\
                             ('nKs', '{}int16'.format(maxCalcs)),('ns', '{}int8'.format(maxCalcs))])
    

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
#local data
iplot = -1
for ipath, path in enumerate(paths): 
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
#         print 'struct',struct
        os.chdir(struct)
        if collateMeshMat:
            for calc in calcs:
                os.chdir(calc)
                if os.path.exists('cell_IBZmesh.m'):
                    meshPlots.write('\t(* {} *)\n'.format(calc))
                    lines = readfile('cell_IBZmesh.m')
                    lines.append( '\n\n')
                    meshPlots.writelines(lines)
                os.chdir('../')          
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
        ns = [] #the base n of the run run
        nDone = 0
        if useSym:
            try:
                nops,IBZvolcut,IBZvol = readSym(calcs[0])
            except:
                sys.exit('Stopping. readSym failed. Set useSym to False')
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
            errs = abs(energies-eref)*1000 + 1e-4 #now in meV 
            data[iplot]['ID'] = '{} {}'.format(struct,tag)
            nAtoms = getNatoms('{}/POSCAR'.format(calc))
            data[iplot]['nAtoms'] = nAtoms
            if useSym:
                data[iplot]['nops'] = nops
                data[iplot]['IBZvolcut'] = IBZvolcut
            data[iplot]['nDone'] = nDone
            data[iplot]['eners'][:nDone] = energies
            data[iplot]['errs'][:nDone] = errs
            data[iplot]['nKs'][:nDone] = nKs
            data[iplot]['ns'][:nDone] = ns
            data[iplot]['color'] = color
            method = path.split('_')[-1].split('/')[0]
            data[iplot]['method'] = method
        os.chdir(path) 
# os.chdir(extpath)


if not extpaths is None:
    for ipath,extpath in enumerate(extpaths):
        os.chdir(extpath)
        method = extpath.split('/')[0]
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
        if collateMeshMat:
            meshPlots = open('IBZmeshPlots','w')
            allMeshesLocal = readfile('{}/IBZmeshPlots'.format(paths[0]))                      
        atomdirs = sorted([d for d in os.listdir(extpath) if os.path.isdir(d) and filter in d])# os.chdir(extpath)
        for dir in atomdirs:
            atom = dir.split('_')[0]
            os.chdir(dir)
            structs = []
            for item in os.listdir(os.getcwd()): 
                if os.path.isdir(item):
                    if item[0].isdigit():
                        os.system('mv {} {}_{}'.format(item,atom,item.split('_')[0])) #rename to format 'Cu_11' instead of 11_atom
            for item in os.listdir(os.getcwd()): 
                if os.path.isdir(item) and (filter2 == None or filter2 in struct):
                    structs.append(item)
            for struct in structs:
                if collateMeshMat:
                    #read sym operators
                    [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = readposcar('POSCAR',dir)
                #         create_poscar('POSCAR',descriptor, scale, latticevecs, natoms, postype, positions, path) #just to remove the scale problem
                    os.chdir(dir)
                    totatoms = sum(natoms)
                    atype = 1
                    aTypes = []
                    for natom in natoms:
                        for i in range(natom):
                            aTypes.append(atype)
                        atype += 1
                    aTypes = array(aTypes)
                    [symopsList, fracsList] = get_spaceGroup(transpose(latticevecs),aTypes,transpose(positions),1e-3,postype.lower()[0] == 'd')
                    self.nops = len(symopsList)
                    self.symops = zeros((3,3,self.nops),dtype = float)
                    for iop in range(len(symopsList)):
                        self.symops[:,:,iop] = trimSmall(array(symopsList[iop]))
                    bounds = [[],[]]
                    blines = readfile('bounds')
                    for line in blines:
                        bounds[0].append([float(str) for str in lines.split()[:3]])
                        bounds[1].append(lines.split()[4]) 
                    for i, line in allMeshesLocal:
                        if struct in line:
                            ilineStruct = i
                            break
                os.chdir(struct)
                iplot += 1
                energies = []
                nKs = []
                ns = [] #the base n of the run run
                nDone = 0
                if useSym:
#                     try:
                    nops,IBZvolcut,nAtoms = copyData(structfile,data)
#                     except:
#                         sys.exit('Stopping. copyData failed. Set useSym to False')
                for calc in calcs:  
                    if electronicConvergeFinish(calc):
                        ener = getEnergy(calc) #in energy/atom
                        if not areEqual(ener,0,1e-5):
                            nDone +=1
                            energies.append(ener)
#                             if 'vc' in path:
                            nK = getNkIBZ(calc,'KPOINTS')
#                             else:
#                                 nK = getNkIBZ(calc,'IBZKPT')
                            if nK > maxNk: maxNk = nK
                            nKs.append(nK)
                            ns.append(int(calc.split('_')[-1]))
                    if collateMeshMat:
                        for j,line in allMeshesLocal[ilineStruct:]:
                            if calc in line.split():
                                IBZfacets = allMeshesLocal[ilineStruct+j+1]
                                meshPlots.write(IBZfacets)
                                break
                        klines = readfile('KPOINTS')
                        ravg = (det(reciplatt)/float(nK))**(1/3.0)
                        eps = ravg/200
                        meshPoints = zeros((nK,3),dtype = float)
                        for i,line in enumerate(klines[3:3+nK]):
                            meshPointDirect0 = [float(string) for string in line.split()[:3]]
                            meshPointDirect1 = meshPointDirect0
                            for j in meshPointDirect1: #get into first BZ
                                while meshPointDirect1[j] > 0.5: 
                                    meshPointDirect1[j] -= 1.0
                                while meshPointDirect1[j] < -0.5:
                                    meshPointDirect1[j] += 1.0
                            meshPoint0 = cartFromDirect(meshPointDirect1)
                            for iop in range(nops):
                                meshPoint = dot(self.symops[:,:,iop],meshPoint0)
                                if isOutside(meshPoint,bounds,eps):
                                    continue
                                else:
                                    break
                            else:
                                sys.exit('Symmetry operations did not bring kpoint {} {} {} into the IBZ'.format(meshPointDirect0[0],meshPointDirect0[1],meshPointDirect0[2]))   
                            meshPoints[i,:] = meshPoint  

                            def facetsMeshMathFile(self,cell,tag,color):
                                '''Output for Mathematica graphics drawing BZ facets and spheres at each mesh point'''
                                strOut = ''
                                strOut = self.facetsMathToStr(strOut,cell,'s','True','Red'); 
                                strOut += ';\n p=Graphics3D[{'
                                cell.mesh = list(trimSmall(array(cell.mesh)))
                                for ipoint,point in enumerate(cell.mesh):
                                    if color != None:
                                        strOut += '{},'.format(color)
                                    strOut += 'Opacity[0.7],Sphere[{' + '{:12.8f},{:12.8f},{:12.8f}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(self.rpacking)
                                    if ipoint < len(cell.mesh) -1:
                                        strOut += ','
                                strOut += '}];\nShow[s,p]'
                                writefile(strOut,'cell_{}.m'.format(tag))   

                            
                                
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
                    errs = abs(energies-eref)*1000 + 1e-4 #now in meV 
                    data[iplot]['ID'] = '{} {}'.format(struct,tag)
                    nAtoms = getNatoms('{}/POSCAR'.format(calc))
                    data[iplot]['nAtoms'] = nAtoms
                    if useSym:
                        data[iplot]['nops'] = nops
                        data[iplot]['IBZvolcut'] = IBZvolcut
                    data[iplot]['nDone'] = nDone
                    data[iplot]['eners'][:nDone] = energies
                    data[iplot]['errs'][:nDone] = errs
                    data[iplot]['nKs'][:nDone] = nKs
                    data[iplot]['ns'][:nDone] = ns
                    data[iplot]['color'] = color
                    method = path.split('_')[-1].split('/')[0]
                    data[iplot]['method'] = method

                os.chdir('../')
            os.chdir(extpath)
        if collateMeshMat:
            meshPlots.close()


nplots = iplot+1 

lines = [' ID , nKIBZ , ener , err, nAtoms, nops,IBZcut\n']  
for iplot in range(nplots):
    n = data[iplot]['nDone']
    for icalc in range(n):#data[iplot]['eners'][:n].tolist()
        lines.append('{}_n{},{},{:15.12f},{:15.12f},{},{},{}\n'.format(data[iplot]['ID'],\
          data[iplot]['ns'][icalc], data[iplot]['nKs'][icalc],\
          data[iplot]['eners'][icalc],data[iplot]['errs'][icalc],\
          data[iplot]['nAtoms'],data[iplot]['nops'],data[iplot]['IBZvolcut']))
writefile(lines,'{}/summary.csv'.format(summaryPath)) 

#plots
if filter[0] == '_':filter = '' #labels can't begin with _
# plotTypes = ['linear','loglog'] #loglinear
# print 'plot only loglog'
plotTypes = ['loglog'] #loglinear
# plotTypes = [] 
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
    methods2 = []
    for iplot in range(nplots):
        labelStr = None
        n = data[iplot]['nDone']
        if coloring == 'method':  
            method = data[iplot]['method'] 
            data[iplot]['color'] = colorsList[methods.index(method)] 
            if method != oldmethod and method not in methods2:
                if doLabel: labelStr = '{} {}'.format(filter,data[iplot]['method'])
                plotData(data[iplot],n,plotType,filter,doLegend,labelStr)
                oldmethod = method;labelStr = None
                methods2.append(method)
            else:
                plotData(data[iplot],n,plotType,filter,doLegend,labelStr)
        elif coloring == 'indiv': 
            if doLabel: labelStr = '{} {}'.format(filter,data[iplot]['ID'])
            plotData(data[iplot],n,plotType,filter,doLegend,labelStr)
#Method averaging
if coloring == 'method':
    print 'Averaging, plotting method errors'
    nbins = int(10*ceil(log10(maxNk)))# 10 bins per decade
    nKstart = 4
    nKbins = array([(nKstart*(10.0**(1/10.0))**i) for i in range(nbins)])
    print 'Nkbins',nKbins
    fig = figure()
    ax1 = fig.add_subplot(111)
    xlabel('N k-points (smoothed by factor {})'.format(int(smoothFactor)))
    ylabel('Error (meV)') 
    methodCostsLogs = []
    for im,method in enumerate(methods):
        methnKmax = 0 
        binCounts = zeros(nbins,dtype = int32)
        binErrsLog = zeros(nbins,dtype = float)
        binErrsWorst = zeros((nbins,nplots),dtype = float)
        costLogs = zeros(nbins,dtype = float) # "Costs" relative to Si Monkhorst Pack, which has err = 10^3/nK^3 + 10^-3 meV.         
        for iplot in range(nplots):
            if data[iplot]['method'] == method:
                for icalc in range(data[iplot]['nDone']-1):
                    nK = data[iplot]['nKs'][icalc]
                    if nK>methnKmax: methnKmax = nK
                    if nK>=nKstart:
                        for ibin in range(nbins):
                            if abs(log10(nK/nKbins[ibin])) <= log10(smoothFactor)\
                              and nKbins[ibin]<= maxNk:
                                err = data[iplot]['errs'][icalc]
                                if err > binErrsWorst[ibin,iplot]:
                                    binErrsWorst[ibin,iplot] = err
                                binErrsLog[ibin] += log10(err)
                                costLogs[ibin] += log10(err/(10**3/(nK**3.0)+0.001))
                                binCounts[ibin] += 1
#         binErrsWorst2 = binErrsWorst[where(binErrsWorst>0)[0]]
        binErrsWorstMean = array([mean(binErrsWorst[ibin,:]) for ibin in range(nbins)])
        mask = where(binErrsWorstMean>0)
        binErrsWorst2 = binErrsWorstMean[mask[0]]
        mask = where(binCounts>0)
        binErrsLog2 = binErrsLog[mask[0]]
        binCounts2 = binCounts[mask[0]]
        nKbins2 = nKbins[mask[0]]
        costLogs2 = costLogs[mask[0]]
        nbins2 = len(nKbins2)
        avgErrs = [binErrsLog2[ibin]/binCounts2[ibin] for ibin in range(nbins2)]
        avgErrsLins = [10**avgErrs[ibin] for ibin in range(nbins2)]
        avgCostLogs =  [costLogs2[ibin]/binCounts2[ibin] for ibin in range(nbins2)]
        avgCostLins = [10**avgCostLogs[ibin] for ibin in range(nbins2)]
        methodCostsLogs.append(mean(avgCostLogs))
        loglog(nKbins2,avgErrsLins,label = method,\
              color = colorsList[im], marker = None)
        loglog(nKbins2,avgCostLins,label = None,\
              color = colorsList[im], marker = None,linestyle=':')
        loglog(nKbins2,binErrsWorst2,label = None,\
              color = colorsList[im], marker = None,linestyle='dashed')
        print 'Method',method, 'nKmax',methnKmax, 'avgLogCost', mean(avgCostLogs)
#         print 'nKs',nKbins2
#         print 'errs',binErrsLog2
#         print 'counts',binCounts2
        print 'avgErrs',avgErrs
    legend(loc='lower left',prop={'size':12});
    fig.savefig('{}/methodErrs'.format(summaryPath))

if useSym:
    #Method averaging with scaled Nk to account for symmetry and nAtoms advantages
    maxNk = 0
    for iplot in range(nplots):
        nops = data[iplot]['nops'] 
        nAtoms = data[iplot]['nAtoms'] 
        n = data[iplot]['nDone']
        for icalc in range(n):#data[iplot]['eners'][:n].tolist()
            nK = data[iplot]['nKs'][icalc]
            if nK>1:
                data[iplot]['nKs'][icalc] = data[iplot]['nKs'][icalc]*nAtoms*nops
                maxNi = max(data[iplot]['nKs'])
                if maxNi > maxNk: maxNk = maxNi
            
    lines = [' ID , nKIBZ , ener , err, nAtoms, nops,IBZcut\n']  
    for iplot in range(nplots):
        n = data[iplot]['nDone']
        for icalc in range(n):#data[iplot]['eners'][:n].tolist()
            lines.append('{}_n{},{},{:15.12f},{:15.12f},{},{},{}\n'.format(data[iplot]['ID'],\
              data[iplot]['ns'][icalc], data[iplot]['nKs'][icalc],\
              data[iplot]['eners'][icalc],data[iplot]['errs'][icalc],\
              data[iplot]['nAtoms'],data[iplot]['nops'],data[iplot]['IBZvolcut']))
    writefile(lines,'{}/summaryScaled.csv'.format(summaryPath))    
    print '\nAveraging, plotting with nK scaling of symmetry and nAtoms'
    nbins = int(10*ceil(log10(maxNk)))# 10 bins per decade
    nKstart = 4
    nKbins = array([(nKstart*10.0**(1/10.0))**i for i in range(nbins)])
    fig = figure()
    ax1 = fig.add_subplot(111)
    smoothFactor *= 2 #double smoothFactor
    xlabel('N k-points (smoothed by factor {}, scaled by N-atoms and N-sym)'.format(int(smoothFactor)))
    ylabel('Error (meV)')
    methodCostsLogs = []
    for im,method in enumerate(methods):
        methnKmax = 0 
        binCounts = zeros(nbins,dtype = int32)
        binErrsLog = zeros(nbins,dtype = float)
        binErrsWorst = zeros(nbins,dtype = float)
        costLogs = zeros(nbins,dtype = float) # "Costs" relative to Si Monkhorst Pack, which has err = 10^3/nK^3 + 10^-3 meV.         
        for iplot in range(nplots):
            if data[iplot]['method'] == method:
                for icalc in range(data[iplot]['nDone']-1):
                    nK = data[iplot]['nKs'][icalc]
                    if nK>methnKmax: methnKmax = nK
                    if nK>1:
                        for ibin in range(nbins):
                            if abs(log10(nK/nKbins[ibin])) <= log10(smoothFactor)\
                              and nKbins[ibin]<= maxNk:
                                err = data[iplot]['errs'][icalc]
                                if err > binErrsWorst[ibin]:
                                    binErrsWorst[ibin] = err
                                binErrsLog[ibin] += err
                                costLogs[ibin] += log10(err/(10**3/(nK**3.0)+0.001))
                                binCounts[ibin] += 1
        mask = where(binCounts>0)
        binErrsLog2 = binErrsLog[mask[0]]
        binErrsWorst2 = binErrsWorst[mask[0]]
        binCounts2 = binCounts[mask[0]]
        nKbins2 = nKbins[mask[0]]
        costLogs2 = costLogs[mask[0]]
        nbins2 = len(nKbins2)
        avgErrs = [binErrsLog2[ibin]/binCounts2[ibin] for ibin in range(nbins2)]
        avgCostLogs =  [costLogs2[ibin]/binCounts2[ibin] for ibin in range(nbins2)]
        avgCostLins = [10**avgCostLogs[ibin] for ibin in range(nbins2)]
        methodCostsLogs.append(mean(avgCostLogs))
        loglog(nKbins2,avgErrs,label = method,\
              color = colorsList[im], marker = None)
        loglog(nKbins2,avgCostLins,label = None,\
              color = colorsList[im], marker = None,linestyle=':')

        print 'Method',method, 'nKmax',methnKmax, 'avgLogCost', mean(avgCostLogs)
#         print 'nKs',nKbins2
#         print 'errs',binErrsLog2
#         print 'counts',binCounts2
#         print 'avgErrs',avgErrs
    legend(loc='lower left',prop={'size':12});
    fig.savefig('{}/methodErrsScaledNk'.format(summaryPath))
print 'Done'