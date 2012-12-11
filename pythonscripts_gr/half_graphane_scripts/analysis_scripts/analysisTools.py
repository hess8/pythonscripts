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
        print('Testing OSZICAR in: ' + i+'\n')
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
    print folder
    proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    try:
        numions = int(newstring[0].split()[-1])
    except:
        return 100 # can't read distance
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
    return mindistance

def writeCCDistances(checkedList):
    lastfolder = os.getcwd()
    '''write distances of C-C expansions to file'''
    ccdistfile = open('ccdistances','w')
    for ielement,ipath in enumerate(checkedList):
        newdist = getCCDistance(ipath)  
        ratio = newdist #/1.46724  #for graphene
        ccdistfile.write(str(ratio) +'\n')
    ccdistfile.close()
    os.chdir(lastfolder)  
        
def getCCDistance(folder):
    os.chdir(folder)
    proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    try:
        numions = int(newstring[0].split()[-1])
    except:
        return 100 # can't read distance
    proc3 = subprocess.Popen(['grep','-i','A3','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc3.communicate()
    repeat = math.fabs(float(newstring[0].split()[-1].split(")")[0]))
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    try:
        line = proc2.communicate()[-2].split('\n')[-2].split(':')[0]
    except:
        return 100 #can't read distance
    outcar = open('OUTCAR','r')
    text = outcar.readlines()
    outcar.close()
    carbon1=[float(text[int(line)+1].split()[0]),float(text[int(line)+1].split()[1]),float(text[int(line)+1].split()[2])]
    carbon2=[float(text[int(line)+2].split()[0]),float(text[int(line)+2].split()[1]),float(text[int(line)+2].split()[2])]
    distance1 = math.sqrt(math.pow(carbon1[2]-carbon2[2],2))
    carbon1[2] = carbon1[2]-repeat
    distance2 = math.sqrt(math.pow(carbon1[2]-carbon2[2],2))
    carbon1[2] = carbon1[2] + repeat
    carbon2[2] = carbon2[2] - repeat
    distance3 = math.sqrt(math.pow(carbon1[2]-carbon2[2],2))
    mindistance =  min(distance1,distance2,distance3)
    return mindistance

             
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
   