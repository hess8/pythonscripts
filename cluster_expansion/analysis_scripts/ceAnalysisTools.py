import os, string, subprocess, math, numpy as np, matplotlib as p

#class fits:
#    def __init__(self,fitsDir,paths,runArray):
#        self.FitsDir = fitsDir
#        self.Paths = paths #those with completed fits
#        self.N = len(paths)
#        self.RunArray = runArray 
#        self.n2 = 0 
#        self.n3 = 0 
#        self.n4 = 0 
#        self.n5 = 0 
#        self.n6 = 0       
   
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
            
def fillRunArray(checkedList, varsList):
    '''Multidimensional array with results of runs:  runArray[nstruc, n2body,growvar] '''
    nComplete = 0
    structureslist = [int(i) for i in varsList[0]]
    clusterlist = [int(i) for i in varsList[1]]
    growlist = [round(float(i),2) for i in varsList[2]]
    print structureslist
    print clusterlist
    print growlist
    dim_nstruc = len(structureslist)  #these need to match run dimensions
    dim_n2body = len(clusterlist)
    dim_growvar = len(growlist)   
    runArray = np.zeros((dim_nstruc,dim_n2body,dim_growvar,5), dtype=float)
    lowestErr = 1000.0
    paths = []
    print checkedList
    for path in checkedList:
        [nstruc, nfits, n2, growvar] = getValues(path)
#        print [nstruc, n2, growvar]
        os.chdir(path)
        [avgErr,stdevErr,L1,L0] = [0,0,0,0]
#        try:
        resultsfile = open('results.out','r')
        results = resultsfile.readlines()[1:] #remove header
        if len(results) == nfits:
            nComplete += 1
            paths.append(path)       
            try:
                os.system('date > complete.txt')
                [avgErr,stdevErr,L1,L0] = resultsOut(results) #over the nfits cases
                print [avgErr,stdevErr,L1,L0]
                if avgErr < lowestErr:
                    lowestErr = avgErr
                    bestfit = [path,[nstruc, n2, growvar],[avgErr,stdevErr,L1,L0]]
                i1 = structureslist.index(nstruc)
                i2 = clusterlist.index(n2)
                i3 = growlist.index(growvar)
    #                    print i1,i2,i3
                runArray[i1,i2,i3,0]=avgErr
                runArray[i1,i2,i3,1]=stdevErr
                runArray[i1,i2,i3,2]=L1
                runArray[i1,i2,i3,3]=L0                     
            except: 
                print 'failed to analyze %s' % [nstruc, n2, growvar]                      
        else:
            print 'results.out length is %s in [nstruc, n2, growvar]: %s' % (len(results),[nstruc, n2,growvar])
#            print avgErr, stdevErr            
#        except:
        print 'no results.out in %s' % [nstruc, n2, growvar]
#    print runArray[2,3,:,1]   
#    plotArray(runArray)     
    print 'number of incomplete jobs', len(checkedList)-nComplete
    print 'Best fit', bestfit
    print 'runArray done'
    return [runArray, paths, bestfit]

def recordFits(runArray,paths,bestfit,maxclust):
    '''read J.1.out: order,avg distance,vertices,J's
    test clusters, identify new ones, assign index to each, 
    record cluster indices for fit'''
    bestpath = bestfit[0]
    bestindices = bestfit[1]
    [bestorders,bestJs,bestverts,bestdists] = readJ1(bestpath,maxclust)
#    print  'best',  [bestorders,bestJs,bestverts,bestdists]
#    mag =        
    clIndex = 0
    for path in paths:
        print path
        dot = 0.0
        [orders,js,verts,dists] = readJ1(bestpath,maxclust)
        for i,jEn in enumerate(js): #check "inner product" between this fit and the best fit
#            print 'i,jEn',i,jEn           
            for ncl in range(orders[i]):
                [match, ibest] = checkCluster(i,orders[i],verts[i,:,:],dists[i],bestorders,bestverts, bestdists)
#                if match:
#                    dot = dot + jEn*bestJs[ibest]          
        
def checkCluster(i,orderi,vertsi,disti,bestorders,bestverts, bestdists):   
    '''determines whether this cluster matches one of the best fit clusters, and returns the index in the best fit'''
    if i == 0:
        match = True #all have constant term
        ib = i
        print 'cluster %i matches bestfit cluster %i' % (i,ib)
        return [match,ib]
    match = False
    ib = 1
#    print 'i',i
#    for ib,orderb in enumerate(bestorders): \
    while ib < len(bestorders) and match == False:
#        print ib
        orderb =  bestorders[ib]
        if orderb == orderi:
#            print orderb, orderi
            for dist in bestdists: #check all avg distances
#                print disti , dist
                if abs(disti - dist) > 1e-6: # not a match 
                    continue #to next distance
                else: #possible match                    
                    vectmatch=np.zeros((orderb), dtype = int)
                    for ivert in range(orderb):
                        bvert = bestverts[ib,ivert,:3]
#                        print bvert, vertsi[ivert]
                        if np.array_equal(bvert,vertsi[ivert]) or np.array_equal(bvert,-vertsi[ivert]): #any new vertex means no match
                            vectmatch[ivert] = 1
#                        print vectmatch
                    if np.sum(vectmatch) == orderb:
                        match = True
                        imatch = ib
                        break  #don't continue to other distances                 
        ib += 1 
    print 'Cluster %s in current fit, order %s, matches cluster %s  in bestfit, %s' % (str(i), str(orderi),str(imatch), str(match))
    return [match,ib]
                    
def readJ1(path,maxclust):
    '''Reads all clusters that are included in the fit (number is L0 norm).  Each cluster has its order and its vertices
    stored in an array'''
#    runArray = np.zeros((L0max,L0max,L0max,4), dtype=float)   
    maxclust = 220
    orders = np.zeros((maxclust),dtype=int)
    verts = np.zeros((maxclust,6,3), dtype=float) #vertices
    dists = np.zeros((maxclust), dtype=float) #vertices
    js = np.zeros((maxclust), dtype=float)
    j1file = open(path+'J.1.out','r')
    lines = j1file.readlines()
    j1file.close()
    index = 0
    verts[index,0,0]=0.0 #1-body cluster
    verts[index,0,1]=0.0
    verts[index,0,2]=0.0 
    orders[0]=1
    print orders[:20]
    js[0] = lines[4].split()[1]
    for i,line in enumerate(lines):
    #    print i, line
        if 'Number of vertices' in line:
            index += 1        
#            print i, lines[i+1]
            order = int(lines[i+1])
            orders[index] = order
            print orders[:20]
            js[index] = float(lines[i-1].split()[0])
            dists[index] = float(lines[i+3])
#            print "order", order
            for j in range(order):
                [x,y,z] = lines[i+7+j].split()[:3] #x, y, z coordinates
                verts[index,j,0]=float(x)
                verts[index,j,1]=float(y)
                verts[index,j,2]=float(z) 
#    print  orders
#    print js
#    print verts
#    print dists                 
    return [orders,js,verts,dists]

def plotArray(x,y,matrix1,plotfile1,title1,xlabel1,ylabel1,plotmax):
    '''plots colored matrix for 2 d array'''
#    from __future__ import division
    from matplotlib.patches import Patch
    from pylab import *
    print plotfile1
#    print x, y
#    x=np.append(x,x[-1])#duplicate last value in extra slot, so plot will show all rows/columns 
#    y=np.append(y,y[-1])
#    print x,y  
    X,Y = meshgrid(x, y)
    Z = matrix1
    fig = figure()
    pcolor(X, Y, Z, cmap=cm.hot, vmax = plotmax)
    xlim((x.min(),x.max()))
    title(title1)
    xlabel(xlabel1)
    ylabel(ylabel1)
    colorbar()
    show()
    fig.savefig(plotfile1)

#pylab.ylabel('voltage (mV)')
#pylab.title('About as simple as it gets, folks')    
    
            
def readList(listname):
    file1 = open(listname,'r')
    list1 = nstrip(file1.readlines())       
    file1.close()
    return list1

def resultsOut(list1):
    '''Gets avg, stdev, L1, L0 norms from columns of results.out'''
    err = 0.0
    stdev = 0.0
    L1 = 0.0
    L0 = 0      
    for line in list1:
        err += float(line.split()[3])
        L1 += float(line.split()[6])
        L0 += float(line.split()[7])               
    err = err/len(list1)
    L1=L1/len(list1)
    L0=L0/len(list1)
    for line in list1:
        stdev += (float(line.split()[3])-err)**2
    stdev = math.sqrt(stdev)/len(list1)
    return [err, stdev,L1,L0]

def getValues(path): 
    '''gets each tag's value from the path'''
    for segment in path.split('/'):
        testsplit = segment.split('_')
        if 'nfitstruc' in testsplit:
            nstruc = int(testsplit[1])
        if 'nfits' in testsplit:
            nfits = int(testsplit[1])                
        if 'n2body' in testsplit:
            n2 = int(testsplit[1])                        
        if 'n3body' in testsplit:
            n3 = int(testsplit[1]) 
        if 'n4body' in testsplit:
            n4 = int(testsplit[1]) 
        if 'n5body' in testsplit:
            n5 = int(testsplit[1]) 
        if 'n6body' in testsplit:
            n6 = int(testsplit[1])
        if 'grow' in testsplit:  
            growVar = round(float(testsplit[1]),2)
            n3 = int(n2 * growVar)
            n4 = int(n3 * growVar)                                                                                       
            n5 = int(n4 * growVar)
            n6 = int(n5 * growVar)
    return [nstruc, nfits, n2, growVar] 
           
def writeEnergiesOszicar(checkedList):           
    lastfolder = os.getcwd()
    enerfile = open('energies','w')
    for i in checkedList:
        os.chdir(i)
        try:
            oszicar = open(i+'OSZICAR','r')          
            energy = oszicar.readlines()[-1].split()[2]
            oszicar.close()
        except:
    		energy = '0'
        enerfile.write(energy + '\n') #energy in last line
    enerfile.close()
    os.chdir(lastfolder) 
        
def writeDistances(checkedList,structure):
    '''write distances of adatoms to file'''
    lastfolder = os.getcwd()
    distfile = open('distances','w')
    for ielement,ipath in enumerate(checkedList):
        distfile.write(str(getDistance(ipath,structure)) +'\n')
    distfile.close()
    os.chdir(lastfolder)

def readResultsOut(folder,structure): 
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

def writeCCDistances(checkedList):
    '''write distances of C-C expansions to file'''
    lastfolder = os.getcwd()
    ccdistfile = open('ccdistances','w')
    diffzfile = open('diffz','w')
    for ielement,ipath in enumerate(checkedList):
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

def nstripfloat(list):
#    '''Strips off /n and converts to float'''
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
    
def elConvergeCheck(folder,NSW):
    """Tests electronic convergence is done by whether the electronic step is less than NELM."""
    try:
        value = getElSteps(folder)
        return value < NSW #True/False
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
    
def writeElSteps(checkedList):    
    '''writes number of steps to output file for each folder'''
    stepsfile = open('elsteps','w')
    for ielement,path in enumerate(checkedList):
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
    
    
def writeSteps(checkedList):    
    '''writes number of steps to output file for each folder'''
    stepsfile = open('steps','w')
    for ielement,path in enumerate(checkedList):
        stepsfile.write(str(getSteps(path))+'\n')
    stepsfile.close()
    
    
def writeFinish(checkedList): 
    '''Writes Y or N depending on vasp finishing, for runs other than relaxation'''
    finishfile = open('finish','w')
    for ielement,path in enumerate(checkedList):
        #get element name     
        if FinishCheck(path):
            finishfile.write('Y' +'\n')
        else:
            finishfile.write('N' +'\n')
    finishfile.close()
   
def writeConverge(checkedList): 
    '''Writes Y or N depending on convergence AND vasp finishing'''
    convergefile = open('converge','w')
    #get NSW, the max ionic steps allowed in the run.  Using first directory in checkedList
    proc = subprocess.Popen(['grep','-i','NSW',checkedList[0]+'/INCAR'],stdout=subprocess.PIPE)
    NSW = int(proc.communicate()[0].split('=')[-1])
    for ielement,path in enumerate(checkedList):
        #get element name
#        element = getElement('adatom_',path)
#        print element      
        if convergeCheck(path,NSW) and FinishCheck(path):
            convergefile.write('Y' +'\n')
        else:
            convergefile.write('N' +'\n')
    convergefile.close()
    
def writeElConverge(checkedList): 
    '''Writes Y or N depending on convergence AND vasp finishing'''
    elconvergefile = open('elconverge','w')
    #get NELM, the max electronic steps allowed in the run.  Using first directory in checkedList
    proc = subprocess.Popen(['grep','-i','NELM',checkedList[0]+'/INCAR'],stdout=subprocess.PIPE)
    result =  proc.communicate()[0]
    NELM = int(result.split('=')[1].split()[0])
    file1 = open('nelm','w')
    file1.write(str(NELM))
    file1.close()
    for ielement,path in enumerate(checkedList):  
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



   