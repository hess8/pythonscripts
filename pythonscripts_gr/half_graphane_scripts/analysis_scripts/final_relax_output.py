#print the name of files to analyze
#Specify Directory to use
mainDir = '/bluehome/bch/vasprun/graphene.structures/half_graphane/'

#Specify the subdir
subdir = 'final_relax'

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
checkedList=checkFolders(toCheckList,checkedList,run)

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

################# spreadsheet #################
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


outfile = open('isolated_atoms.csv','w')
outfile.write('Element,Calculated Energy,Distance, CC Distance\n')
for i in range(len(elements)):
    linei = elements[i]+','+energies[i]+','+distances[i]+','+ccdistances[i]+'\n'
    outfile.write(linei)
outfile.close()

print 'Done'


