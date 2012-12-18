#print the name of files to analyze
#Specify Directory to use
mainDir = '/bluehome/bch/vasprun/graphene.structures/half_graphane/'

#Specify the subdir
subdir = 'isolated'

dir = mainDir + subdir + '/'
#Specify the name of the type of run
runName = 'relaxation' 

import os,subprocess,math,time 
import numpy as np 
from analysisTools import addToList, checkFolders, writeEnergiesOszicar, getElement, writeElements, nstrip

run = runName

#            
print('\nInitializing...\n')
print('Finding Directories to do %s on\n' % run)
print('Searching all Directories in ' + dir+'\n')
toCheckList = [] #these inits must be here because of recursion in functions using them
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


################# spreadsheet #################
#Open data files
os.chdir(dir)
outfile = open('isolated_atoms.csv','w')

file = open('elements','r')  
elements = nstrip(file.readlines())
file.close()
print elements

file = open('energies','r')
energies = nstrip(file.readlines())
file.close()

outfile.write('Element,Calculated Energy\n')
for i in range(len(elements)):
    linei = elements[i]+','+energies[i]+'\n'
    outfile.write(linei)
outfile.close()

print "done"


