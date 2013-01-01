#### Bug: doesn't submit jobs on the same run as it creates folders.  Must run twice if starting 
# with an empty directory.

#print the name of files to analyze
#Specify Directory to use
mainDir = '/bluehome/bch/vasprun/graphene.structures/h.half_graphane/'
#mainDir = "/bluehome/bch/vasprun/graphene.structures/ds_diam_like/"


#Specify the subdir
subdir = 'final_relax'
#subdir = 'test2'
dir = mainDir + subdir + '/'
#Specify the name of the type of run
runName = 'relax' 

#Specify Poscar variables
poscarVariables = {
}

#Specify KPoints variables
kpointVariables = {
}

#Specify Incar variables
incarVariables = {
'@ISPIN':
	['2'],
'@IBRION':
	['2'],
'@ISIF':
	['4']
}

import os,sys,subprocess,math,time 
import numpy as np 
sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/analysis_scripts')
from analysisTools import addToList, checkFolders, writeEnergiesOszicar,  \
    writeElements, nstrip, writeDistances, writeCCDistances, writeConverge, \
    FinishCheck


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

#write out energies from all elements
writeEnergiesOszicar(checkedList)

#Write distance of adatoms
writeDistances(checkedList)

#Write C-C expansion
writeCCDistances(checkedList)

#Check convergence
writeConverge(checkedList) 

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

file = open('converge','r')
converged = nstrip(file.readlines())
file.close()

outfile = open('final_relax.csv','w')
outfile.write('Element,Calculated Energy,Distance,CC Distance,CC Diffz,Converged\n')
for i in range(len(elements)):
    linei = elements[i]+','+energies[i]+','+distances[i]+','+ccdistances[i]+','+diffz[i]+','+converged[i]+'\n'
    outfile.write(linei)
outfile.close()

print 'Done with summary'

################# spreadsheet #################
#Bring in data from other runs
os.chdir(mainDir)
file = open('/bluehome/bch/vasprun/graphene.structures/half_graphane/isolated/energies','r')
mainDir 
isolenergies = nstrip(file.readlines())
file.close()

try:
	file = open('initial_relax/stretch','r')
	strenergies = nstrip(file.readlines())
except:
	strenergies = ['']*len(elements)
file.close()

eIsolatedH = -1.115
eIsolatedC = -1.3179
energyGraphane = -25.63
bindEnergyGraphane = energyGraphane - 2*eIsolatedH - 2* eIsolatedC

binde = [0]*len(elements)
for i in range(len(elements)):
	try:
		#for half graphane:
#		benergy_system = float(energies[i]) - float(isolenergies[i]) - eIsolatedH - 2*eIsolatedC 
#		binde[i] = benergy_system - bindEnergyGraphane #relative to graphane

        # for h.half_graphane:
#		benergy_system = float(energies[i]) - float(isolenergies[i]) - 2*eIsolatedH - 2*eIsolatedC 
		binde[i] = float(energies[i]) - energyGraphane - float(isolenergies[i]) #energy to bring adatoms from far away to beneath H layer of graphane
		#		binde[i] = float(energies[i]) - 2*float(isolenergies[i])-2* eIsolatedC 
#		binde[i] = benergy_system - bindEnergyGraphane #relative to graphane
		#For binding E of adatom relative to 
		#for double sided (no H) diamond_like, not relative
#		binde[i] = float(energies[i]) - 2*float(isolenergies[i])-2* eIsolatedC 		
	except:
		binde[i] = 100 #can't read energy
#	if elements[i] == 'Ti':
#		print float(energies[i]) , float(isolenergies[i]), binde[i]
outfile = open('analysis.csv','w')
outfile.write('Element,BE.adatom.underh,Calc Energy,Distance,CC Diffz,CC expans %,Stretch energy,Converged\n')
# write spreadsheet

for i in range(len(elements)):
    try:
        ccexpand = str(round(100*(float(ccdistances[i])/1.53391 - 1),1)) #compare to graphane 
    except:
        ccexpand = 'null'
    if converged[i] =='Y':
         linei = elements[i]+','+str(binde[i])+','+energies[i]+','+distances[i]+','+diffz[i]+','+ccexpand+','+strenergies[i]+','+converged[i]+'\n'
    else:
         linei = elements[i]+'*,'+str(binde[i])+','+energies[i]+'*,'+distances[i]+'*,'+diffz[i]+','+ccexpand+'*,'+strenergies[i]+'*,'+converged[i]+'\n'        
    outfile.write(linei)
outfile.close()

print 'Done with analysis'


