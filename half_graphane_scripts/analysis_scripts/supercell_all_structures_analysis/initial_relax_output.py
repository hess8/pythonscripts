import os,subprocess,math,time 
from numpy import array, zeros, binary_repr,log2
from analysisTools import FinishCheck,allFoldersList,energyOszicar,nadatoms


maindir = '/fslhome/bch/cluster_expansion/hexagonal/2x2adatoms/'
run = 'relax'
enerIsolatedAd=	-.45365181E+01 # for tungsten (W) only
nsites = 8 

#standard energies:
eIsolatedH = -1.115
eIsolatedC = -1.3179
energyGraphane = -25.63
nPrimCells = nsites/2
bindEnergyGraphane = nPrimCells * (energyGraphane - 2*eIsolatedH - 2* eIsolatedC)

dirlist = allFoldersList(maindir,run)
nstruct = len(dirlist)
#print dirlist
ener = zeros(nstruct)
structs = []
distances = zeros(nstruct)
binaryatoms = zeros(nstruct)
nadtom = zeros(nstruct)
for i,dir in enumerate(dirlist):
    os.chdir(dir)
    struct = dir.split('/')[-3]
    structs.append(struct)
    ener[i] = energyOszicar()
    binaryatoms[i] = binary_repr(int(struct.replace('struct','')))
    nAd[i] = nadatoms(str(binaryatoms[i]))
    nH = nsites-nadtom[i]
    benergy_system = ener[i] -  nAd[i] * enerIsolatedAd - nH*eIsolatedH - nsites*eIsolatedC 
    binde[i] = benergy_system - bindEnergyGraphane #relative to graphane 
    BEString = 'BE.vs.graphane'
#    print i, structs[i], ener[i]

os.chdir(maindir)
outfile = open('2x2all_initial.csv','w')
outfile.write('Structure,binary,%s,Calculated Energy\n' % BEString)

#for i in range(len(elements)):

#outfile.close()
for i,dir in enumerate(dirlist):
    linei = structs[i]+','+ str(binaryatoms[i])+','+ binde[i] +','+ str(ener[i])+','+ nadtom[i]+ ','+ BindEnergy[i]+'\n'
    outfile.write(linei)

outfile.close()




#
##Find distance of adatom for minimum energy
#os.chdir(dir)
#resultsfile = open('allenergies','r')
#results = resultsfile.readlines()
#ndist = len(poscarVariables['@distance'])
#nelements = len(results)/ndist
#energies = [0]*ndist
#for ielement in range(nelements):
##    print('ielement %d \n' % ielement)
#    elementstart = ielement*ndist
#    #get element name
#    path = checkedList[elementstart] # just the first folder in list
#    prefix = 'adatom_'
#    index1 = path.index(prefix) 
#    index2 = path.index('/',index1)
#    element = path[index1+len(prefix):index2]
#    print ('Element ' + element)
#    #get energies
#    elementresults = results[elementstart:elementstart+ndist]
##    if element == 'Ti':
##        print 'Ti energies', elementresults
##    print 'elementresults',elementresults
#    for i in range(len(elementresults)):
#        try:
#            energies[i]=float(elementresults[i])
#        except:
#            energies[i]=0.0 #if text, etc
#    enerfar = energies[ndist-1] #the last energy (sorted checkedlist) is the farthest away, most likely to be finished
#    for i in range(len(energies)):
#        if abs(energies[i]-enerfar) > 100: #throw away outliers
#            energies[i] = 0
#    minindex = np.argmin(energies)
##    if element == 'Ti':
##        print 'Ti energies', energies, minindex
##    print ('best energy %f' % energies[minindex])
#    
#    enerstretch = enerfar - energies[minindex]
#    
#    print ('best index %d' % minindex)
#    #get distance from OUTCAR
#    bestfolder = checkedList[elementstart+minindex]
#    bestpathfile.write(bestfolder+'\n')
#    print ('getDistance from %s' % bestfolder )
#    print getDistance(bestfolder)
#    elemfile.write(element +'\n')
#    distfile.write(np.str(getDistance(bestfolder)) +'\n')
#    enerfile.write(np.str(energies[minindex]) +'\n')
#    strchfile.write(np.str(enerstretch) +'\n')
#    
#resultsfile.close()
#elemfile.close()
#enerfile.close()
#distfile.close()
#strchfile.close()
#bestpathfile.close() 
#
############# Spreadsheet 
#import os, sys, string
#sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/')
#from ScriptTools import nstrip
#
#os.chdir(dir)
#outfile = open('half_graphane_initial.csv','w')
#
#file = open('elements','r')  
#elements = nstrip(file.readlines())
#file.close()
#print elements
#
#file = open('energies','r')
#energies = nstrip(file.readlines())
#file.close()
#
#file = open('stretch','r')
#stretch = nstrip(file.readlines())
#file.close()
#
#file = open('distances','r')
#distances = nstrip(file.readlines())
#file.close()
#
#outfile.write('Element,Calculated Energy,Stretch Energy,Distance\n')
#for i in range(len(elements)):
#    linei = elements[i]+','+energies[i]+','+stretch[i]+','+distances[i]+'\n'
#    outfile.write(linei)
#outfile.close()
#
#print 'Done'

