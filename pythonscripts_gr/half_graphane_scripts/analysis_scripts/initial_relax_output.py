#print the name of files to analyze
#Specify Directory to use
mainDir = '/bluehome/bch/vasprun/graphene.structures/half_graphane/'

#Specify the subdir
subdir = 'initial_relax'

dir = mainDir + '/' + subdir + '/'
#Specify the name of the type of run
runName = 'relaxation' 

#Specify Poscar variables
poscarVariables = {
'@distance':
	[8,6,4,3,2,1]  
}

import os,subprocess,math,time 
import numpy as np 


run = runName
newRun = 'DOS'
newRunFile = 'DOSCAR' #Will check to see if higher level is already done, then skip it

global toCheckList
global checkedList
global toRunList

def addToList(folder):
    files = os.listdir(folder)
#    print files
    for path in files:
#        print path
        if os.path.isdir(folder+path+'/'):
            toCheckList.append(folder+path+'/')
            addToList(folder+path+'/')
#            print path+'/'
#
def checkFolders():
    for path in toCheckList:
#        print('CHECK NEXT LINE')
#        print(path.split('/')[-2])
        if path.split('/')[-2] == run:
            checkedList.append(path)
#            
def getDistance(folder):
    lastfolder = os.getcwd()
    os.chdir(folder)
    print folder
    proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    numions = int(newstring[0].split()[-1])
    proc3 = subprocess.Popen(['grep','-i','A3','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc3.communicate()
    repeat = math.fabs(float(newstring[0].split()[-1].split(')')[0]))
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    try:
        line = proc2.communicate()[-2].split('\n')[-2].split(':')[0]
    except:
        return 100 #can't read distance
    outcar = open('OUTCAR','r')
    text = outcar.readlines()
    outcar.close()
    carbon1=[float(text[int(line)+1].split()[0]),float(text[int(line)+1].split()[1]),float(text[int(line)+1].split()[2])]
    adatom1=[float(text[int(line)+3].split()[0]),float(text[int(line)+3].split()[1]),float(text[int(line)+3].split()[2])]
    distance1 = math.sqrt(math.pow(carbon1[0]-adatom1[0],2)+math.pow(carbon1[1]-adatom1[1],2)+math.pow(carbon1[2]-adatom1[2],2))
    adatom1[2] = adatom1[2]-repeat
    distance2 = math.sqrt(math.pow(carbon1[0]-adatom1[0],2)+math.pow(carbon1[1]-adatom1[1],2)+math.pow(carbon1[2]-adatom1[2],2))
    adatom1[2] = adatom1[2] + repeat
    carbon1[2] = carbon1[2] - repeat
    distance3 = math.sqrt(math.pow(carbon1[0]-adatom1[0],2)+math.pow(carbon1[1]-adatom1[1],2)+math.pow(carbon1[2]-adatom1[2],2))
    mindistance =  min(distance1,distance2,distance3)
    os.chdir(lastfolder) #restore to last folder
    return mindistance
    

print('\nInitializing...\n')
print('Finding Directories to do %s on\n' % run)
print('Searching all Directories in ' + dir+'\n')
toCheckList = []
checkedList = []
toRunList = []
addToList(dir)

print toCheckList
print('Checking to see which folders contain '+run+'\n')
time.sleep(1)
checkFolders()

print '\nThe following folders are in checkedList:'
for i in checkedList:
    print('checkedList contains : ' + i+'\n')
    
#write out energies from all runs
os.chdir(dir)
file = open('allenergies','w')
for i in checkedList:
    os.chdir(i)
    print('Testing OSZICAR in: ' + i+'\n')
    oszicar = open(i+'OSZICAR','r')
    energy = oszicar.readlines()[-1].split()[2]
    file.write(energy + '\n') #energy in last 
    oszicar.close()
file.close()
os.chdir(dir)
#Open output files
elemfile = open('elements','w')
enerfile = open('energies','w')
distfile = open('distances','w')
strchfile = open('stretch','w')
bestpathfile = open('bestpath','w')


#Find distance of adatom for minimum energy
os.chdir(dir)
resultsfile = open('allenergies','r')
results = resultsfile.readlines()
ndist = len(poscarVariables['@distance'])
nelements = len(results)/ndist
energies = [0]*nelements
for ielement in range(nelements):
    print('ielement %d \n' % ielement)
    elementstart = ielement*ndist
    #get element name
    path = checkedList[elementstart] # just the first folder in list
    prefix = 'adatom_'
    index1 = path.index(prefix) 
    index2 = path.index('/',index1)
    element = path[index1+len(prefix):index2]
    print ('Element ' + element)
    
    #get energies
    elementresults = results[elementstart:elementstart+ndist]
    print 'elementresults',elementresults
    for i in range(len(elementresults)):
        try:
            energies[i]=float(elementresults[i])
        except:
            energies[i]=0.0 #if text, etc
    enerfar = energies[0] #the first energy is the farthest away, most likely to be finished
    for i in range(len(energies)):
        if abs(energies[i]-enerfar) > 100: #throw away outliers
            energies[i] = 0
    minindex = np.argmin(energies)
    print ('best energy %f' % energies[minindex])
    
    enerstretch = enerfar - energies[minindex]
    
    print ('best index %d' % minindex)
    #get distance from OUTCAR
    bestfolder = checkedList[elementstart+minindex]
    bestpathfile.write(bestfolder+'\n')
    print ('getDistance from %s' % bestfolder )
    print getDistance(bestfolder)
    elemfile.write(element +'\n')
    distfile.write(np.str(getDistance(bestfolder)) +'\n')
    enerfile.write(np.str(energies[minindex]) +'\n')
    strchfile.write(np.str(enerstretch) +'\n')
    
resultsfile.close()
elemfile.close()
enerfile.close()
distfile.close()
strchfile.close()
bestpathfile.close() 

print 'Done'
