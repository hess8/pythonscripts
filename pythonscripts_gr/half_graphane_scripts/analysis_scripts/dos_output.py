#Specify Directory to use
#mainDir = '/bluehome/bch/vasprun/graphene.structures/h.half_graphane2.1/'
#mainDir = "/bluehome/bch/vasprun/graphene.structures/ds_diam_like/"
#mainDir = "/bluehome/bch/vasprun/graphene.structures/transmet.half_graphane/h.half_graphane2.1/"
mainDir = "/bluehome/bch/vasprun/graphene.structures/transmet.half_graphane/half_graphane/"

isolatedDir = '/bluehome/bch/vasprun/graphene.structures/transmet.half_graphane/isolated/electron_relax/'

#Specify the subdir
subdir = 'dos'
#subdir = 'test1'
dir = mainDir + subdir + '/'
#Specify the name of the type of run
runName = 'dos' 

#get type of structure
lastdir = mainDir.split('/')[-2]
if 'h.' in lastdir:
    structure = 'h.'
elif 'diam' in lastdir:
    structure = 'diam' #doublesided diamondlike
else: #half_graphane
    structure = 'half_gr'

import os,sys,subprocess,math,time
import matplotlib.pyplot as plt 
import numpy as np 
sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/analysis_scripts')
from analysisTools import addToList, checkFolders, writeEnergiesOszicar,  \
    writeElements, nstrip, writeDistances, writeCCDistances, \
    finishCheck, writeElSteps, writeElConverge, getElement

run = runName

print('\nInitializing...\n')
print('Finding Directories to do %s on\n' % run)
print('Searching all Directories in ' + dir+'\n')
toCheckList = []
checkedList = []
 
#find directories
os.chdir(dir)
toCheckList = addToList(dir,toCheckList)

print('Checking to see which folders contain '+run+'\n')
time.sleep(1)
checkedList=sorted(checkFolders(toCheckList,checkedList,run))

print '\nThe following folders are in checkedList:'
for i in checkedList:
    print('checkedList contains : ' + i+'\n')
   
#write out elements list    
writeElements(checkedList)
print 'elements', os.getcwd() 
#write out energies from all elements
writeEnergiesOszicar(checkedList)
print 'energies', os.getcwd() 
#Write distance of adatoms
writeDistances(checkedList, structure)
print 'distances', os.getcwd() 
#Write C-C expansion
writeCCDistances(checkedList)
print 'CCdist', os.getcwd() 
#Check convergence
writeElConverge(checkedList) 
print 'ElConverge', os.getcwd() 
#number of relaxation steps
writeElSteps(checkedList)
print 'writeElSteps', os.getcwd() 
###### plot DOS #####
sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/analysis_scripts/plotting/')
print 'Plotting DOS'
from plotTools import vasputil_dosplot, collate_plots, plotArray
for path in checkedList:
    try:
        vasputil_dosplot([], ["DOSCAR"], path) #options, args, dir
        plt.close()
    except:
        print 'Fail:'+ path
print 'Collating plots'
os.chdir(mainDir)
tablelist = []
tableorder = ['Ti','V','Cr','Mn','Fe','Co','Ni',
              'Zr','Nb_pv','Mo','Tc','Ru','Rh','Pd',
              'Hf','Ta','W','Re','Os','Ir','Pt']
for elem in tableorder:
    for path in checkedList:
        if elem in path:
            tablelist.append(path)
                
collate_plots(tablelist,'dos.png','dos.png')

################# summary spreadsheet #################
#Open data files

os.chdir(dir)

file = open('elements','r')  
elements = nstrip(file.readlines())
file.close()
print elements

file = open('energies','r')
energies = nstrip(file.readlines())
file.close()

file = open('distances','r')
distances = nstrip(file.readlines())
file.close()

file = open('ccdistances','r')
ccdistances = nstrip(file.readlines())
file.close()

file = open('diffz','r')
diffz = nstrip(file.readlines())
file.close()

file = open('elconverge','r')
converged = nstrip(file.readlines())
file.close()

file = open('elsteps','r')
steps = nstrip(file.readlines())
file.close()

outfile = open('dos.csv','w')
outfile.write('Element,Calculated Energy,Distance,CC Distance,CC Diffz,Converged\n')
for i in range(len(elements)):
    linei = elements[i]+','+energies[i]+','+distances[i]+','+ccdistances[i]+','+diffz[i]+','+converged[i]+'\n'
    outfile.write(linei)
outfile.close()

print 'Done with summary'


################# spreadsheet #################
#Bring in data from other runs
os.chdir(mainDir)
file = open(isolatedDir+'energies','r')
isolenergies = nstrip(file.readlines())
file.close()

file = open(isolatedDir+'convdiff','r')
isolconvdiff = nstrip(file.readlines())
file.close()

file = open(isolatedDir+'finish','r')
isolatedFinish = nstrip(file.readlines())
file.close()

file = open(isolatedDir+'elsteps','r')
isolatedElSteps = nstrip(file.readlines())
file.close()

file = open(isolatedDir+'nelm','r')
NELM = nstrip(file.readlines())[0]
file.close()

try:
    file = open('stretch/energies','r')
    strenergies = nstrip(file.readlines())
except:
    strenergies = ['']*len(elements)
file.close()

eIsolatedH = -1.115
eIsolatedC = -1.3179
energyGraphane = -25.63
bindEnergyGraphane = energyGraphane - 2*eIsolatedH - 2* eIsolatedC
BEstretch = [0]*len(elements)
binde = [0]*len(elements)
for i in range(len(elements)):  
    BEstretch[i] = float(strenergies[i])-float(energies[i])
for i in range(len(elements)):  
    try:
        # for h.half_graphane:
        if structure == 'h.':
            binde[i] = float(energies[i]) - energyGraphane - float(isolenergies[i]) #energy to bring adatoms from far away to beneath H layer of graphane
            BEString = 'BE.adatom.underh'
#not used:        benergy_system = float(energies[i]) - float(isolenergies[i]) - 2*eIsolatedH - 2*eIsolatedC 

        #for double sided (no H) diamond_like, not relative
        elif structure == 'diam':
            binde[i] = float(energies[i]) - 2*float(isolenergies[i])-2* eIsolatedC     
            BEString = 'BE absolute'
        #for half graphane:
        elif structure == 'half_gr': 
            benergy_system = float(energies[i]) - float(isolenergies[i]) - eIsolatedH - 2*eIsolatedC 
            binde[i] = benergy_system - bindEnergyGraphane #relative to graphane 
            BEString = 'BE.vs.graphane'
        else:
            print "UNKNOWN STRUCTURE"
            BEString = 'unknown structure'
    except:
        binde[i] = 100 #can't read energy or structure not right
#    if elements[i] == 'Ti':
#        print float(energies[i]) , float(isolenergies[i]), binde[i]
outfile = open('dos_analysis.csv','w')
outfile.write('Element,'+BEString+',Calc Energy,Isol atom,IsolConvDiff,Distance,CC Diffz,CC expans %,BEstretch,Stretch energy,Converged,Steps\n')
# write spreadsheet

for i in range(len(elements)):
#    print i, elements[i], isolenergies[i], isolatedFinish[i] 
    try:
        ccexpand = str(round(100*(float(ccdistances[i])/1.53391 - 1),1)) #compare to graphane 
    except:
        ccexpand = 'null'
    if converged[i] =='Y' and isolatedFinish[i] == "Y" and float(isolatedElSteps[i]) < NELM:
        linei = elements[i]+','+str(binde[i])+','+energies[i]+','+isolenergies[i]+','+isolconvdiff[i]+','+distances[i]+','+diffz[i]+','+ccexpand+','+str(BEstretch[i])+','+strenergies[i]+','+converged[i]+','+steps[i]+'\n'
    elif isolatedFinish[i] == "N" or float(isolatedElSteps[i]) == NELM:
        print i, elements[i],isolatedFinish[i], float(isolatedElSteps[i]) 
        linei = elements[i]+'*,'+str(binde[i])+','+energies[i]+'*,'+'not done'+','+isolconvdiff[i]+','+distances[i]+'*,'+diffz[i]+','+ccexpand+'*,'+str(BEstretch[i])+','+strenergies[i]+'*,'+converged[i]+'*,'+steps[i]+'\n'        
    else:
        linei = elements[i]+'*,'+str(binde[i])+','+energies[i]+'*,'+isolenergies[i]+'*,'+isolconvdiff[i]+','+distances[i]+'*,'+diffz[i]+','+ccexpand+'*,'+str(BEstretch[i])+','+strenergies[i]+'*,'+converged[i]+'*,'+steps[i]+'\n'        
    print i,linei
    outfile.write(linei)
outfile.close()

################# data matrices plotting #################
# Put in matrices of periodic table shape, and correct order
nrows = 3
ncols = 7
enerMat = np.zeros((nrows+1,ncols+1),dtype = float)#need extra so will show all data
distMat = np.zeros((nrows+1,ncols+1),dtype = float)
diffzMat = np.zeros((nrows+1,ncols+1),dtype = float)
ccexpMat = np.zeros((nrows+1,ncols+1),dtype = float)

for i in range(len(tableorder)):
    row = np.floor(i/(ncols))
    col = np.mod(i,ncols)
    print i, [row, col]
    for j in range(len(elements)):
        if tableorder[i] == elements[j]:
            print elements[j]
            if binde[j] < 0:
                enerMat[row,col] = -binde[j]
            else:
                enerMat[row,col] = 0.0
            distMat[row,col] = distances[j]
            diffzMat[row,col] = diffz[j]
            ccexpMat[row,col] = round(100*(float(ccdistances[j])/1.53391 - 1),1) #compare to graphane
#    #switch row order
#    for i in range(len(enerMat[:,1])):

#for i in range(nrows):
#    print enerMat[i,:]
        
#plotting matrices

os.chdir(mainDir)
x = np.asarray(range(ncols+1))
y = -np.asarray(range(nrows+1))
#print x,y
matrix = enerMat
plotmax = np.amax(matrix)
plotfile = 'bindenergy_%s' % structure
title1 = 'Binding energy for %s' % structure
xlabel1 = ''
ylabel1 = ''
plotArray(x,y,matrix,plotfile,title1,xlabel1,ylabel1,plotmax)
                        

print 'Done with analysis'