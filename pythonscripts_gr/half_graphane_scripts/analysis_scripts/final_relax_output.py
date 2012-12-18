#### Bug: doesn't submit jobs on the same run as it creates folders.  Must run twice if starting 
# with an empty directory.

#print the name of files to analyze
#Specify Directory to use
mainDir = '/bluehome/bch/vasprun/graphene.structures/half_graphane/'

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

import os,subprocess,math,time 
import numpy as np 
from analysisTools import addToList, checkFolders, writeEnergiesOszicar,  \
    writeElements, nstrip, writeDistances, writeCCDistances


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


outfile = open('final_relax.csv','w')
outfile.write('Element,Calculated Energy,Distance,CC Distance\n')
for i in range(len(elements)):
    linei = elements[i]+','+energies[i]+','+distances[i]+','+ccdistances[i]+'\n'
    outfile.write(linei)
outfile.close()

print 'Done with summary'

################# spreadsheet #################
#Bring in data from other runs
os.chdir(mainDir)
file = open('isolated/energies','r')
isolenergies = nstrip(file.readlines())
file.close()

file = open('initial_relax/stretch','r')
strenergies = nstrip(file.readlines())
file.close()

eIsolatedH = -1.115
eIsolatedC = -1.3179
energyGraphane = -25.63
bindEnergyGraphane = energyGraphane - 2*eIsolatedH - 2* eIsolatedC

binde = [0]*len(elements)
for i in range(len(elements)):
	try:
		binde[i] = float(energies[i]) - float(isolenergies[i])-eIsolatedH - 2* eIsolatedC -bindEnergyGraphane 
	except:
		binde[i] = 100 #can't read energy
	if elements[i] == 'Ti':
		print float(energies[i]) , float(isolenergies[i]), binde[i]
outfile = open('half_graphane_analysis_run_again_after_finalRelax2.csv','w')
outfile.write('Element,Binding Energy,Calculated Energy,Distance,CC expansion, Stretch energy\n')

# write spreadsheet
for i in range(len(elements)):
	#compare ccdistances to graphane
    linei = elements[i]+','+str(binde[i])+','+energies[i]+','+distances[i]+','+str(float(ccdistances[i])/1.467)+','+strenergies[i]+'\n'
    outfile.write(linei)
outfile.close()

ccratio = str(float(ccdistances[i])/1.467) #compare to graphane

print 'Done with analysis'


