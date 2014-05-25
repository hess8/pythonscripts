import os, string, subprocess, math
from numpy import zeros

def readfile(filepath):
    file1 = open(filepath,'r')
    lines = file1.readlines()
    file1.close()
    return lines

def writefile(lines,filepath):
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

def addToList(folder,toCheckList):
    files = os.listdir(folder)
    for path in files:
        if os.path.isdir(folder+path+'/'):
            toCheckList.append(folder+path+'/')
            addToList(folder+path+'/',toCheckList)
    return toCheckList

def checkFolders(toCheckList,dirslist,run):
    dirslist = []
    for path in toCheckList:
        if path.split('/')[-2] == run:
            dirslist.append(path)
    return dirslist            
            

def writedirnames(list):    
    namesfile = open('names','w')
    for name in list:
        namesfile.write(name +'\n')
    namesfile.close()

def writeEnergiesOszicar(list):           
    lastfolder = os.getcwd()
    enerfile = open('energies','w')
    for i in list:
        print i
#        os.chdir(i)
#        print os.getcwd()
#        print os.listdir(os.getcwd())
        try:
            oszicar = open(i+'/OSZICAR','r')
            energy = oszicar.readlines()[-1].split()[2]
            oszicar.close()
        except:
    		energy = '0'
        print energy
        enerfile.write(energy + '\n') #energy in last line
        os.chdir(lastfolder)
    enerfile.close()
    os.chdir(lastfolder)
     
def enerparts(list):
    '''Finds the 9 terms contributing to the total energy in vasp, as in:
         Free energy of the ion-electron system (eV)
          ---------------------------------------------------
          alpha Z        PSCENC =        -0.43400971
          Ewald energy   TEWEN  =      -150.19868154
          -1/2 Hartree   DENC   =        -0.23270923
          -exchange  EXHF       =         0.00000000
          -V(xc)+E(xc)   XCENC  =       -52.00285543
          PAW double counting   =        71.87483513       -8.01524693
          entropy T*S    EENTRO =        -0.00011188
          eigenvalues    EBANDS =        24.48126738
          atomic energy  EATOM  =       107.07526000'''
 
    enerparts = zeros((len(list),9),dtype = float)
    for ifolder, folder in enumerate(list):
        try:
            outcar = open(folder + '/OUTCAR','r')
            text = outcar.readlines()
            outcar.close()           
            for i in range(len(text)): #find last line matching text
                if 'Free energy of the ion-electron system' in text[i]:
                    nline = i
            for i in range(9):
                enerparts[ifolder,i] = float(text[nline+2+i].split('=')[-1].split()[0]) #to handle PAW line with two entries
        except:
            'go on' # can't read them
    return enerparts
    
    
def writeNkIBZ(list):     
    '''number of K points extracted from IBZKPT'''
    lastfolder = os.getcwd()
    file = open('NkIBZ','w')
    for i in list:
        print i
#        os.chdir(i)
#        print os.getcwd()
#        print os.listdir(os.getcwd())
        try:
            ibzkpt = open(i+'/IBZKPT','r')
#            print ibzkpt.readlines()[1].split()
            nk = ibzkpt.readlines()[1].split()[0]
            ibzkpt.close()
        except:
            nk = '0'
        print nk
        file.write(nk + '\n')  
        os.chdir(lastfolder)
    file.close()
    os.chdir(lastfolder) 
    
def writeNk(list):     
    '''number of K points intended,from comment in KPOINTS'''
    lastfolder = os.getcwd()
    file = open('Nk','w')
    for i in list:
        print i
#        os.chdir(i)
#        print os.getcwd()
#        print os.listdir(os.getcwd())
        try:
            ibzkpt = open(i+'/KPOINTS','r')
#            print ibzkpt.readlines()[1].split()
            nk = ibzkpt.readlines()[0].split()[0]
            ibzkpt.close()
        except:
            nk = '0'
        print nk
        file.write(nk + '\n')  
        os.chdir(lastfolder)
    file.close()
    os.chdir(lastfolder)
        
def writeDistances(dirslist,structure):
    '''write distances of adatoms to file'''
    lastfolder = os.getcwd()
    distfile = open('distances','w')
    for ielement,ipath in enumerate(dirslist):
        distfile.write(str(getDistance(ipath,structure)) +'\n')
    distfile.close()
    os.chdir(lastfolder)

def getDistance(folder,structure): 
    '''gets adatom-carbon distance'''
    os.chdir(folder)
    try:
        outcar = open('OUTCAR','r')
        text = outcar.readlines()
        proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        numions = int(newstring[0].split()[-1])
    except:
        return 100 # can't read distance
    proc3 = subprocess.Popen(['grep','-n','lattice vectors','OUTCAR'],stdout=subprocess.PIPE)
    nline = proc3.communicate()[-2].split('\n')[-2].split(':')[0] #returns one line after grep
    try:
        repvector=[float(text[int(nline)+2].split()[0]),float(text[int(nline)+2].split()[1]),float(text[int(nline)+2].split()[2])]
        repeat = repvector[2]
    except:
        repeat = 100.0
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    try:
        nline = proc2.communicate()[-2].split('\n')[-2].split(':')[0] 
    except:
        return 100 #can't read distance
    outcar.close()
    if structure == 'h.':
        adatomline = 5           
    elif structure == 'diam':
        adatomline = 4     #actually 2 adatoms here, could be 3 or 4
    elif structure == 'half_gr': #half_graphane
        adatomline = 4
    else:
        print "UNKNOWN STRUCTURE"
    carbon1line = int(nline)+1   
    carbon2line = int(nline)+2
    adatomline = int(nline)+adatomline
    carbon1=[float(text[carbon1line].split()[0]),float(text[carbon1line].split()[1]),float(text[carbon1line].split()[2])]
    carbon2=[float(text[carbon2line].split()[0]),float(text[carbon2line].split()[1]),float(text[carbon2line].split()[2])]
    adatom1=[float(text[adatomline].split()[0]),float(text[adatomline].split()[1]),float(text[adatomline].split()[2])]
    distancemin = 100
    for i in [-1,1]:
        for j in [-1,1]:
            for carbon in [carbon1, carbon2]:
                carbontry = [carbon[0], carbon[1], carbon[2] + i*repeat]
                adatomtry = [adatom1[0], adatom1[1], adatom1[2] + j*repeat]
                distancenew = distance(carbontry,adatomtry)
                if distancenew < distancemin:
                    distancemin = distancenew              
    return distancemin


def distance(vec1, vec2):
	return math.sqrt(math.pow(vec1[0]-vec2[0],2)+math.pow(vec1[1]-vec2[1],2)+math.pow(vec1[2]-vec2[2],2))

def writeCCDistances(dirslist):
    '''write distances of C-C expansions to file'''
    lastfolder = os.getcwd()
    ccdistfile = open('ccdistances','w')
    diffzfile = open('diffz','w')
    for ielement,ipath in enumerate(dirslist):
        newdist = getCCDistance(ipath)
        ccdistfile.write(str(newdist[0]) +'\n')
        diffzfile.write(str(newdist[1]) +'\n')
    ccdistfile.close()
    diffzfile.close()
    os.chdir(lastfolder)  
        
def getCCDistance(folder):
    os.chdir(folder)
    try:
        outcar = open('OUTCAR','r')
        text = outcar.readlines()
        proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        numions = int(newstring[0].split()[-1])
    except:
        return [100,100] # can't read distance
    proc3 = subprocess.Popen(['grep','-n','lattice vectors','OUTCAR'],stdout=subprocess.PIPE)
    nline = proc3.communicate()[-2].split('\n')[-2].split(':')[0] #returns one line after grep
    try:
        repvector=[float(text[int(nline)+2].split()[0]),float(text[int(nline)+2].split()[1]),float(text[int(nline)+2].split()[2])]
        repeat = repvector[2]
    except:
        repeat = 100.0
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    try:
        nline = proc2.communicate()[-2].split('\n')[-2].split(':')[0]
    except:
        return [100,100] #can't read distance
    outcar.close()
    carbon1=[float(text[int(nline)+1].split()[0]),float(text[int(nline)+1].split()[1]),float(text[int(nline)+1].split()[2])]
    carbon2=[float(text[int(nline)+2].split()[0]),float(text[int(nline)+2].split()[1]),float(text[int(nline)+2].split()[2])]
    distance1 = distance(carbon1,carbon2)
    diffz1 = abs(carbon1[2] - carbon2[2])
    carbon1[2] = carbon1[2]-repeat
    distance2 = distance(carbon1,carbon2)
    diffz2 =  abs(carbon1[2] - carbon2[2])
    carbon1[2] = carbon1[2] + repeat
    carbon2[2] = carbon2[2] - repeat
    distance3 = distance(carbon1,carbon2)
    diffz3 =  abs(carbon1[2] - carbon2[2])
    return [min(distance1, distance2, distance3), min(diffz1,diffz2,diffz3)]
             
def getElement(prefix,path):
    '''Prefix is e.g. adatom_'''   
    index1 = path.index(prefix) 
    index2 = path.index('/',index1)
    element = path[index1+len(prefix):index2]
    return element
    
def writeElements(dirslist):    
    elemfile = open('elements','w')
    for ielement,ipath in enumerate(dirslist):
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
    
def convergeCheck(folder,NSW):
    """Tests whether force convergence is done by whether the last line of Oszicar is less than NSW."""
    try:
        value = getSteps(folder)
        return value < NSW #True/False
    except:
        return False #True/False
    
def elConvergeCheck(folder,NELM):
    """Tests electronic convergence is done by whether the electronic step is less than NELM."""
    try:
        value = getElSteps(folder)
        return value < NELM #True/False
    except:
        return False #True/False

def getElSteps(folder):
    '''number of electronic steps for isolated runs'''
    lastfolder = os.getcwd()
    os.chdir(folder)
    try:
        oszicar = open('OSZICAR','r') 
        laststep = oszicar.readlines()[-2].split()[1] # Next to last line, 2nd entry
        oszicar.close()
        os.chdir(lastfolder) 
        value = int(laststep)
        return value         
    except:
        os.chdir(lastfolder)         
        return 9999
    
def writeElSteps(dirslist):    
    '''writes number of steps to output file for each folder'''
    stepsfile = open('elsteps','w')
    for ielement,path in enumerate(dirslist):
        stepsfile.write(str(getElSteps(path))+'\n')
    stepsfile.close()

def getSteps(folder):
    '''number of steps in relaxation, as an integer'''
    lastfolder = os.getcwd()
    os.chdir(folder)
    if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
        os.chdir(lastfolder) 
        return -9999
    oszicar = open('OSZICAR','r')
    laststep = oszicar.readlines()[-1].split()[0]
    oszicar.close()
    os.chdir(lastfolder)  
    try:
        value = int(laststep)
        return value
    except:
        return 9999
    
    
def writeSteps(dirslist):    
    '''writes number of steps to output file for each folder'''
    stepsfile = open('steps','w')
    for ielement,path in enumerate(dirslist):
        stepsfile.write(str(getSteps(path))+'\n')
    stepsfile.close()
    
    
def writeFinish(dirslist): 
    '''Writes Y or N depending on vasp finishing, for runs other than relaxation'''
    finishfile = open('finish','w')
    for ielement,path in enumerate(dirslist):
        #get element name     
        if FinishCheck(path):
            finishfile.write('Y' +'\n')
        else:
            finishfile.write('N' +'\n')
    finishfile.close()
   
def writeConverge(dirslist): 
    '''Writes Y or N depending on convergence AND vasp finishing'''
    convergefile = open('converge','w')
    #get NSW, the max ionic steps allowed in the run.  Using first directory in dirslist
    proc = subprocess.Popen(['grep','-i','NSW',dirslist[0]+'/INCAR'],stdout=subprocess.PIPE)
    NSW = int(proc.communicate()[0].split('=')[-1])
    for ielement,path in enumerate(dirslist):
        #get element name
#        element = getElement('adatom_',path)
#        print element      
        if convergeCheck(path,NSW) and FinishCheck(path):
            convergefile.write('Y' +'\n')
        else:
            convergefile.write('N' +'\n')
    convergefile.close()
    
def writeElConverge(dirslist): 
    '''Writes Y or N depending on convergence AND vasp finishing'''
    elconvergefile = open('elconverge','w')
    #get NELM, the max electronic steps allowed in the run.  Using first directory in dirslist
    proc = subprocess.Popen(['grep','-i','NELM',dirslist[0]+'/INCAR'],stdout=subprocess.PIPE)
    result =  proc.communicate()[0]
    NELM = int(result.split('=')[1].split()[0])
    file1 = open('nelm','w')
    file1.write(str(NELM))
    file1.close()
    for ielement,path in enumerate(dirslist):  
        if elConvergeCheck(path,NELM) and FinishCheck(path):
            elconvergefile.write('Y' +'\n')
        else:
            elconvergefile.write('N' +'\n')
    elconvergefile.close()

def FinishCheck(folder):
#        """Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR."""
    lastfolder = os.getcwd()
    os.chdir(folder)
    proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    os.chdir(lastfolder)    
    return newstring[0].find('Voluntary') > -1 #True/False

def cpuTime(folder):
#    proc = subprocess.Popen(['tail', 'OUTCAR','|','grep','User'],shell=True,stdout=subprocess.PIPE)
    lastfolder = os.getcwd()
    os.chdir(folder)
    proc = subprocess.Popen(['tail OUTCAR | grep User'],shell=True,stdout=subprocess.PIPE)
    result =  proc.communicate()
    os.chdir(lastfolder) 
    try:
        time =  float(result[0][-10:-3].strip()) #last few characters
    except:
        time = 0
    print time
    return time

def writeCPUtime(list):    
    '''writes number of steps to output file for each folder'''
    file = open('cputime','w')
    for path in list:
        print path
        file.write(str(cpuTime(path))+'\n')
    file.close()
    
def getdata(list,name):
    'reads files that have simply a single line. Returned as string.  convert to float, int if needed'
    data = []
    for dir in list:
        data.append(readfile(dir+'/'+name)[0])
    return data

def getms(list):
    '''finds the m factor out of the path name which is in the form f45_8.  8 is the m factor'''
    ms = []
    for dir in list:
        ms.append(dir.split('_')[-1].strip())
    return ms