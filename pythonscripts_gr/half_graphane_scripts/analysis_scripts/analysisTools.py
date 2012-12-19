import os, string, subprocess, math

def addToList(folder,toCheckList):
    files = os.listdir(folder)
    for path in files:
        if os.path.isdir(folder+path+'/'):
            toCheckList.append(folder+path+'/')
            addToList(folder+path+'/',toCheckList)
    return toCheckList

def checkFolders(toCheckList,checkedList,run):
    checkedList = []
    for path in toCheckList:
        if path.split('/')[-2] == run:
            checkedList.append(path)
    return checkedList            
            
def writeEnergiesOszicar(checkedList):           
    lastfolder = os.getcwd()
    enerfile = open('energies','w')
    for i in checkedList:
        os.chdir(i)
#        print('Testing OSZICAR in: ' + i+'\n')
        oszicar = open(i+'OSZICAR','r')
        try:
        	energy = oszicar.readlines()[-1].split()[2]
        except:
    		energy = '0'
        enerfile.write(energy + '\n') #energy in last line
        oszicar.close()
    enerfile.close()
    os.chdir(lastfolder) 
        
def writeDistances(checkedList):
    '''write distances of adatoms to file'''
    lastfolder = os.getcwd()
    distfile = open('distances','w')
    for ielement,ipath in enumerate(checkedList):
        distfile.write(str(getDistance(ipath)) +'\n')
    distfile.close()
    os.chdir(lastfolder)

def getDistance(folder): 

    os.chdir(folder)
    outcar = open('OUTCAR','r')
    text = outcar.readlines()
    proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    try:
        numions = int(newstring[0].split()[-1])
    except:
        return 100 # can't read distance
    proc3 = subprocess.Popen(['grep','-n','lattice vectors','OUTCAR'],stdout=subprocess.PIPE)
    nline = proc3.communicate()[-2].split('\n')[-2].split(':')[0] #returns one line after grep
    repvector=[float(text[int(nline)+2].split()[0]),float(text[int(nline)+2].split()[1]),float(text[int(nline)+2].split()[2])]
    repeat = repvector[2]
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    try:
        nline = proc2.communicate()[-2].split('\n')[-2].split(':')[0] 
    except:
        return 100 #can't read distance
    outcar.close()
    carbon1=[float(text[int(nline)+1].split()[0]),float(text[int(nline)+1].split()[1]),float(text[int(nline)+1].split()[2])]
    adatom1=[float(text[int(nline)+3].split()[0]),float(text[int(nline)+3].split()[1]),float(text[int(nline)+3].split()[2])]
    distance1 = distance(carbon1,adatom1)
    adatom1[2] = adatom1[2]-repeat
    distance2 = distance(carbon1,adatom1)
    adatom1[2] = adatom1[2] + repeat
    carbon1[2] = carbon1[2] - repeat
    distance3 = distance(carbon1,adatom1)
    return min(distance1, distance2, distance3)

def distance(vec1, vec2):
	return math.sqrt(math.pow(vec1[0]-vec2[0],2)+math.pow(vec1[1]-vec2[1],2)+math.pow(vec1[2]-vec2[2],2))

def writeCCDistances(checkedList):
    '''write distances of C-C expansions to file'''
    lastfolder = os.getcwd()
    ccdistfile = open('ccdistances','w')
    for ielement,ipath in enumerate(checkedList):
        newdist = getCCDistance(ipath)  
        ratio = newdist #/1.46724  #for graphene
        ccdistfile.write(str(ratio) +'\n')
    ccdistfile.close()
    os.chdir(lastfolder)  
        
def getCCDistance(folder):
    os.chdir(folder)
    outcar = open('OUTCAR','r')
    text = outcar.readlines()
    proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    try:
        numions = int(newstring[0].split()[-1])
    except:
        return 100 # can't read distance
    proc3 = subprocess.Popen(['grep','-n','lattice vectors','OUTCAR'],stdout=subprocess.PIPE)
    nline = proc3.communicate()[-2].split('\n')[-2].split(':')[0] #returns one line after grep
    repvector=[float(text[int(nline)+2].split()[0]),float(text[int(nline)+2].split()[1]),float(text[int(nline)+2].split()[2])]
    repeat = repvector[2]
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    try:
        nline = proc2.communicate()[-2].split('\n')[-2].split(':')[0]
    except:
        return 100 #can't read distance
    outcar.close()
    carbon1=[float(text[int(nline)+1].split()[0]),float(text[int(nline)+1].split()[1]),float(text[int(nline)+1].split()[2])]
    carbon2=[float(text[int(nline)+2].split()[0]),float(text[int(nline)+2].split()[1]),float(text[int(nline)+2].split()[2])]
    distance1 = distance(carbon1,carbon2)
    carbon1[2] = carbon1[2]-repeat
    distance2 = distance(carbon1,carbon2)
    carbon1[2] = carbon1[2] + repeat
    carbon2[2] = carbon2[2] - repeat
    distance3 = distance(carbon1,carbon2)
    return min(distance1, distance2, distance3)
             
def getElement(prefix,path):    
    index1 = path.index(prefix) 
    index2 = path.index('/',index1)
    element = path[index1+len(prefix):index2]
    return element
    
def writeElements(checkedList):    
    elemfile = open('elements','w')
    for ielement,ipath in enumerate(checkedList):
        #get element name
        element = getElement('adatom_',ipath)
        elemfile.write(element +'\n')
    elemfile.close()
    
def nstrip(list):
#	'''Strips off /n'''
    import string
    list2 = []
    for string1 in list:   
        string2 = string1.strip("\n")
        list2.append(string2)
    return list2

def xxxconvergeCheck(folder,NSW):
#        """Tests whether force converges is done by whether the last line of Oszicar is less than NSW."""
    lastfolder = os.getcwd()
    os.chdir(folder)
    if not os.path.exists('OSZICAR'):
        return False
    oszicar = open('OSZICAR','r')
    lines = oszicar.readlines
    print lines
#    for i1,j1 in enumerate(lines):
#        print  oszicar.readlines()[i]
    print folder
    lastline = oszicar.readlines()[-1].split
    print  lastline
#    lastline = oszicar.readlines()[-1].split()
    oszicar.close()
    os.chdir(lastfolder)      
    try:
        value = int(lastline[0])
        return int(lastline[0])< NSW #True/False
    except:
        return False #True/False
    
def convergeCheck(folder,NSW):
#        """Tests whether force converges is done by whether the last line of Oszicar is less than NSW."""
    lastfolder = os.getcwd()
    os.chdir(folder)
    if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
        return False
    laststep = oszicar.readlines()[-1].split()[0]
    oszicar.close()
    os.chdir(lastfolder)  
    try:
        value = int(laststep[0])
        return int(laststep[0])< NSW #True/False
    except:
        return False #True/False

def writeConverge(checkedList):    
    '''Writes Y or N depending on convergence'''
    convergefile = open('converge','w')
    #get NSW, the max ionic steps allowed in the run.  Using first directory in checkedList
    proc = subprocess.Popen(['grep','-i','NSW',checkedList[0]+'/INCAR'],stdout=subprocess.PIPE)
    NSW = int(proc.communicate()[0].split('=')[-1])
    for element,path in enumerate(checkedList):
        #get element name
        if convergeCheck(path,NSW):
            convergefile.write('Y' +'\n')
        else:
            convergefile.write('N' +'\n')
    convergefile.close()





   