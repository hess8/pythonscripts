#print the name of files to analyze
#Specify Directory to use
mainDir = "/fslhome/bch/cluster_expansion/cluster_size_test"

#Specify the subdir
subdir = 'vary_nclusters_ntrain'

dir = mainDir + subdir + '/'
#Specify the name of the type of run
runName = 'relaxation' 

import os,subprocess,math,time,sys,datetime
import numpy as np 
sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/analysis_scripts')
from analysisTools import addToList, checkFolders, writeEnergiesOszicar,  writeFinish, getElement,\
                            writeElements, nstrip, writeElSteps, writeElConverge

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

#Check number of electronic steps
writeElSteps(checkedList) 

#Check Vasp finish
writeFinish(checkedList) 

#Check Number of electronic steps finish
writeElConverge(checkedList) 

################# spreadsheet #################
#Open data files
os.chdir(dir)
outfile = open('isolated_atoms.csv','w')

file = open('elements','r')  
elements = nstrip(file.readlines())
file.close()

file = open('energies','r')
energies = nstrip(file.readlines())
file.close()

file = open('elsteps','r')
elsteps = nstrip(file.readlines())
file.close()

file = open('elconverge','r')
elconverge = nstrip(file.readlines())
file.close()

file = open('finish','r')
finish = nstrip(file.readlines())
file.close()


if os.path.exists(dir+'oldenergies'):
    file = open('oldenergies','r')
    oldenergies = nstrip(file.readlines())
    file.close()
else:
    oldenergies = ['0.0']*len(elements)       
diffe = ['0.0']*len(elements) 
diffFile = open('convdiff','w')   
outfile.write('Element,Calculated Energy,Diff Ener (conv),Electr Steps, Converge & Finish\n')
for i in range(len(elements)):
#    print elements[i]
#    print energies[i], oldenergies[i]
    try:
        diffe[i] = str(float(energies[i]) - float(oldenergies[i]))
        if float(diffe[i])<1.0e-4:
            diffe[i] = 'done' 
            fileconv = open(checkedList[i]+'converged.txt', 'w') #create a file to indicate it's converged
            fileconv.write('converged')
            fileconv.close()              
    except:
        diffe[i] = '99'
    if finish[i] == 'Y' and elconverge[i] == 'Y':
        converge = 'Y'
    else:
        converge = 'N'
    linei = elements[i]+','+energies[i]+','+diffe[i]+','+elsteps[i]+','+ converge+'\n'
    outfile.write(linei)
    diffFile.write(diffe[i]+'\n')
outfile.close()
diffFile.close()

print 'Done'





