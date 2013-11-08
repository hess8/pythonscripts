import os, math, string, subprocess, sys
from numpy import array, sqrt,zeros,round
from numpy.linalg import norm
from numpy import int as npint
from numpy import int as npfloat

#sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts')
#import poscar #Gus' read poscar class

def allFoldersList(folder,run):
    list1 = []
    list1 = addToList(folder,list1)#lists all dirs and subdirs
    list2 = []
    checkFolders(list1,list2,run) #finds which matches run
    return list2

def addToList(folder,toCheckList):
    files = os.listdir(folder)
    for path in files:
        if os.path.isdir(folder+path+'/'):
            toCheckList.append(folder+path+'/')
            addToList(folder+path+'/',toCheckList)
    return toCheckList  
  
def checkFolders(toCheckList,checkedList,run):
    for path in toCheckList:
        if path.split('/')[-2] == run:
            checkedList.append(path)
    return checkedList            
            
def energyOszicar():           
    try:
        oszicar = open('OSZICAR','r') 
        energy = oszicar.readlines()[-1].split()[2]
        oszicar.close()
    except:
        energy = '0'
    return energy
   
   
def nadatoms(strbinary):
	nad = 0
	for i in range(len(strbinary)):
		nad += int(strbinary[i])
	return nad
	
def FinishCheck():
#        """Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR."""
    proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    return newstring[0].find('Voluntary') > -1 #True/False	

def getSteps():
    '''number of steps in relaxation, as an integer'''
    if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
        return -9999
    oszicar = open('OSZICAR','r')
    laststep = oszicar.readlines()[-1].split()[0]
    try:
        value = int(laststep)
        return value
    except:
        return 9999
    
def getNSW(dir): 
    '''get NSW, the max ionic steps allowed in the run'''   
    lastdir = os.getcwd()
    os.chdir(dir)
    proc = subprocess.Popen(['grep','-i','NSW','INCAR'],stdout=subprocess.PIPE)
    os.chdir(lastdir)
    string = proc.communicate()
    NSW = int(string[0].split('=')[-1])
    return NSW

def directToCartesian(vect,latt):
    '''Convert a vector from a basis in the lattice vectors (direct) to cartesian.
    Also center around the origin by subtracting lattice vectors '''
    for j in range(3):
        if vect[j] > 0.5:
            vect[j] = vect[j] - 1.0 # i
    vect = vect[0]*latt[0,:] + vect[1]*latt[1,:] + vect[2]*latt[2,:] 
    return vect 
     
def getPositionInfo(): 
    ''' Read from CONTCAR '''
    file1 = open('CONTCAR','r'); contcar = file1.readlines(); file1.close()
    #lattice vectors  
    latt = array([[float(string) for string in contcar[2].strip().split()], \
                 [float(string) for string in contcar[3].strip().split()], \
                 [float(string) for string in contcar[4].strip().split()]])
    typelist = contcar[5].strip().split()
    if typelist == ['C','H']:
        [nC, nH] =  [int(string) for string in contcar[6].strip().split()]
        nAd = 0
    elif len(typelist) == 2:  # only two elemennts and not C H, so must be no H
        [nC, nAd] =  [int(string) for string in contcar[6].strip().split()]
        nH = 0
    else:
        [nC, nH, nAd] =  [int(string) for string in contcar[6].strip().split()]
    posC = zeros((nC,3)); posH = zeros((nH,3)); posAd = zeros((nAd,3))
#    print [nC, nH, nAd]
    for i in range(nC):
        posC[i,:] =  [float(string) for string in contcar[8+i].strip().split()]
        posC[i,:] = directToCartesian(posC[i,:],latt)
    for i in range(nH):
        posH[i,:] =  [float(string) for string in contcar[8+nC+i].strip().split()]
        posH[i,:] = directToCartesian(posH[i,:],latt)        
    for i in range(nAd):
        posAd[i,:] =  [float(string) for string in contcar[8+nC+nH+i].strip().split()]
        posAd[i,:] = directToCartesian(posAd[i,:],latt)
    #avg C-C distances
    avgdCC = 0
    avgzC = 0
    count = 0
    for i in range(nC):
        avgzC += abs(posC[i,2])       
        for j in range(nC):
            if i != j:
                if norm(posC[i,:]-posC[j,:])<2.3: #nearest neighbors only (2NN at about >2.5) 
                    count += 1            
                    avgdCC += norm(posC[i,:]-posC[j,:])
#                    print i,j, norm(posC[i,:]-posC[j,:])
    avgdCC = round(avgdCC/count,2); avgzC = round(avgzC/nC,2)
    #Adatom distances
    avgzAd = 0
    avgdAdC = 0
    min_dxyAdC = 100 #initialize for min search
    count = 0
    for i in range(nAd):
        avgzAd += abs(posAd[i,2])       
        for j in range(nC):
            dxyAdC = norm(posAd[i,0:2]-posC[j,0:2])# distance in x-y plane
            if dxyAdC<min_dxyAdC:
                min_dxyAdC = dxyAdC #find minimum xy distance from Adatom to C atom,to see if adatom has moved off top site
            if dxyAdC<0.5*avgdCC:   
                count += 1            
                avgdAdC += norm(posAd[i,:]-posC[j,:])
    if nAd > 0:
        avgdAdC = round(avgdAdC/count,2); avgzAd = round(avgzAd/nAd,2)
    return [avgdCC, avgzC, avgdAdC, avgzAd, round(min_dxyAdC,4)]
