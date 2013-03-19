import os, string, subprocess, math, numpy as np, matplotlib as p
 
   
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
            
def fillRunArray(checkedList, varsList,nHold):
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
    runArray = np.zeros((dim_nstruc,dim_n2body,dim_growvar,6), dtype=float)
    lowestErr = 1000.0
    paths = []
    print checkedList
    for path in checkedList:
        [nstruc, nfits, n2, growvar] = getValues(path)
#        print [nstruc, n2, growvar]
        os.chdir(path)
        [avgErr,stdevErr,L1,L0] = [0,0,0,0]
        try:
            resultsfile = open('results.out','r')
            results = resultsfile.readlines()[1:] #remove header
            if len(results) == nfits:
                nComplete += 1
                paths.append(path)       
                try:
                    os.system('date > complete.txt')
                    [avgErr,stdevErr,L1,L0] = resultsOut(results) #over the nfits cases
#                    print [avgErr,stdevErr,L1,L0]
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
        except:
            print 'no results.out in %s' % [nstruc, n2, growvar]
#    print runArray[2,3,:,1]   
#    plotArray(runArray)  
    print 'number of incomplete jobs', len(checkedList)-nComplete
    print 'Best fit', bestfit
    print 'recording closeness of fits into runArray'
    runArray = recordFits(runArray,paths,bestfit,varsList,np.amax(runArray[:,:,3]),nHold)
    print 'runArray done'
    return [runArray, paths, bestfit]

#def compareFits(paths):
#    '''writes out files for easy comparison'''
    
    
def recordFits(runArray,paths,bestfit,varsList,maxclust,nHold):
    '''read J.1.out: order, avg distance, vertices, J's
    test fits to see how close they are to best fit. '''
    structureslist = [int(i) for i in varsList[0]]
    clusterlist = [int(i) for i in varsList[1]]
    growlist = [round(float(i),2) for i in varsList[2]]    
    bestpath = bestfit[0]
#    print 'structs', structureslist
#    print 'clusters', clusterlist
#    print 'grow',growlist
    bestindices = bestfit[1]
    print 'Best fit parameters'
    [bestorders,bestJs,bestverts,bestdists] = readJ1(bestpath,maxclust)
    bestfitMag = np.sqrt(np.sum(np.square(bestJs)))
#    print  'best',  [bestorders,bestJs,bestverts,bestdists]
    for ipath, path in enumerate(paths):
        print 'Analyzing fit in path', ipath, path
        dot = 0.0
        [orders,js,verts,dists] = readJ1(path,maxclust)
        ordersummary = [(i,orders.tolist().count(i)) for i in range(9)[1:]]
        fitMag = np.sqrt(np.sum(np.square(js)))
#        print 'mag', fitMag
        for i,jEn in enumerate(js): #check "inner product" between this fit and the best fit
            if jEn != 0.0: #using only filled in part of array
                [match, imatch] = checkCluster(i,orders[i],verts[i,:,:],dists[i],bestorders,bestverts, bestdists)
                if match:
#                    print 'i,jEn, imatch, bestJs[imatch]', [[i,jEn],[imatch, bestJs[imatch]]]
                    dot = dot + jEn*bestJs[imatch]   
        closeness = dot/(fitMag*bestfitMag) # essentially cos(theta) = a dot b/(|a||b|)       
#        print 'closeness', closeness
        [nstruc, nfits, n2, growvar] = getValues(path)
        print [nstruc, nfits, n2, growvar]
        os.chdir(path)
        summaryfile = open('fitsummary.out', 'w')
        lines = [path[40:],'closeness to best fit', closeness,nstruc,n2,growvar,'ordersummary',ordersummary,'dists',dists, 'js',js]
        for line in lines:
            summaryfile.write(str(line)+'\n')
        summaryfile.close()
        i1 = structureslist.index(nstruc)
        i2 = clusterlist.index(n2)
        i3 = growlist.index(growvar)
#                    print i1,i2,i3
        runArray[i1,i2,i3,4]=closeness
        
        # prediction errors comparison
        bestPredErrs = readPredErr(bestpath,nHold)
        bestErrMag = np.sqrt(np.sum(np.square(bestPredErrs)))
        runArray[i1,i2,i3,5]=predErrDot(path,bestPredErrs,bestErrMag,nHold)
#        print 'dot', runArray[i1,i2,i3,5]
               
    return runArray

def predErrDot(path,besterrs,bestErrMag,nHold):
    ''' '''
    predErr = readPredErr(path,nHold)
    predErrMag = np.sqrt(np.sum(np.square(predErr)))
    dot = 0.0
    for istr in range(nHold):
        dot += predErr[istr] * besterrs[istr]/predErrMag/bestErrMag
    return dot
#    print getValues(path)

def readPredErr(path,nHold):
    predErr = np.zeros(nHold,dtype=float)
    os.chdir(path)
    file = 'prediction_errors.out'
    fileopen = open(path+file,'r')
    lines = nstrip(fileopen.readlines())[-nHold-1:-1] #get last fit of the Nfit. 
    fileopen.close()
    for i,line in enumerate(lines):
        predErr[i] = float(line.split()[-3]) - float(line.split()[-2])  #[-2], middle column of numbers is VASP data
    return predErr
            
def checkCluster(i,orderi,vertsi,disti,bestorders,bestverts, bestdists):   
    '''determines whether this cluster matches one of the best fit clusters, and returns the index in the best fit'''
    if i == 0:
        match = True #all have constant term
        ib = i
#        print 'cluster %i matches bestfit cluster %i' % (i,ib)
        return [match,ib]
    match = False
    imatch = 0
    ib = 1 #cluster in bestfit to start comparison
#    print 'i',i
#    for ib,orderb in enumerate(bestorders): \
    while ib < len(bestorders) and match == False:
#        print ib
        orderb =  bestorders[ib]
        if orderb == 0:
            break #reached maximum order in best fit.
#        print orderb, orderi        
        if orderb == orderi:
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
#    if match:
#        print 'Cluster %s in fit, order %s, matches cluster %s  in bestfit' % (str(i+1), str(orderi),str(imatch+1))
#    else:
#        print 'No match for cluster %i in fit' % (i+1)
    return [match,imatch]
                    
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
    js[0] = lines[4].split()[1]
    for i,line in enumerate(lines):
    #    print i, line
        if 'Number of vertices' in line:
            index += 1        
#            print i, lines[i+1]
            order = int(lines[i+1])
            orders[index] = order
            js[index] = float(lines[i-1].split()[0])
            dists[index] = float(lines[i+3])
#            print "order", order
            for j in range(order):
                [x,y,z] = lines[i+7+j].split()[:3] #x, y, z coordinates
                verts[index,j,0]=float(x)
                verts[index,j,1]=float(y)
                verts[index,j,2]=float(z) 
#    print  orders
#    print js[:30]
#    print orders[:30]
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

#def readList(listname):
#    file1 = open(listname,'r')
#    list1 = nstrip(file1.readlines())       
#    file1.close()
#    return list1

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
           

def distance(vec1, vec2):
	return math.sqrt(math.pow(vec1[0]-vec2[0],2)+math.pow(vec1[1]-vec2[1],2)+math.pow(vec1[2]-vec2[2],2))

    
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
    
   