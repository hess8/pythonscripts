################## Directories ################## 
mainDir = '/fslhome/bch/cluster_expansion/cluster_size_test/agpt/'
#Specify the type of run
runtype = 'ncl_ntr'
#runtype = 'test'


dir = mainDir + runtype + '/'
#Specify the name of the type of run
runName = 'run_10_15' 

import os,subprocess,math,time,sys,datetime,numpy as np
#import numpy as np 
from ceAnalysisTools import addToList, checkFolders, fillRunArray, readList, plotArray,writeFinish, getElement,\
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

#print '\nThe following folders are in checkedList:'
#for i in checkedList:
#    print('checkedList contains : ' + i+'\n')
print '%s folders will be analyzed' % len(checkedList) 
os.chdir(dir)
structureslist=np.load('structureslist.npy')
clusterlist=np.load('clusterlist.npy')
growlist=np.load('growlist.npy')
#clusterlist=readList('clusterlist.dat')
#growlist=readList('growlist.dat')
varsList = [structureslist,clusterlist,growlist]
runArray = fillRunArray(checkedList, varsList) #also writes complete.txt if results.out has correct Nfits lines. 
print runArray[0,:,:,0]
os.chdir(dir)
################## Plotting ############

plotindex = 0 #[err, stdeverr, L1,L0]
plotmax = np.amax(runArray[:,:,:,plotindex])
x = varsList[2] # choose 0,1,2
y = np.log2(varsList[1])
for i,nstruct in enumerate(structureslist):
    plotfile = 'errstruct%s' % nstruct
    title1 = 'Validation error for %s AgPt structures in training set' % nstruct
    xlabel1 = 'Growth factor for orders above 2-body'
    ylabel1 = 'Log2 of N-pairs'
    plotArray(x,y,runArray[i,:,:,plotindex],plotfile,title1,xlabel1,ylabel1,plotmax)
      
plotindex = 2 #[err, stdeverr, L1,L0]
plotmax = np.amax(runArray[:,:,:,plotindex])
x = varsList[2] # choose 0,1,2
y = np.log2(varsList[1])
for i,nstruct in enumerate(structureslist):
    plotfile = 'l1struct%s' % nstruct
    title1 = 'L1 norm avg for %s AgPt structures in training set' % nstruct
    xlabel1 = 'Growth factor for orders above 2-body'
    ylabel1 = 'Log2 of N-pairs'
    plotArray(x,y,runArray[i,:,:,plotindex],plotfile,title1,xlabel1,ylabel1,plotmax)  

plotindex = 3 #[err, stdeverr, L1,L0]
plotmax = np.amax(runArray[:,:,:,plotindex])
x = varsList[2] # choose 0,1,2
y = np.log2(varsList[1])
for i,nstruct in enumerate(structureslist):
    plotfile = 'l0struct%s' % nstruct
    title1 = 'L0 norm avg for %s AgPt structures in training set' % nstruct
    xlabel1 = 'Growth factor for orders above 2-body'
    ylabel1 = 'Log2 of N-pairs'
    plotArray(x,y,runArray[i,:,:,plotindex],plotfile,title1,xlabel1,ylabel1,plotmax)  
    
plotindex = 0 #[err, stdeverr, L1,L0]
plotmax = np.amax(runArray[:,:,:,plotindex])
x = varsList[2] # choose 0,1,2
y = varsList[0]
for i,npairs in enumerate(clusterlist):
    plotfile = 'npairs%s' % npairs
    title1 = 'Validation error for %s pairs in cluster basis, AgPt ' % npairs
    xlabel1 = 'Growth factor for orders above 2-body'
    ylabel1 = 'Number of structures in training set'
    plotArray(x,y,runArray[:,i,:,plotindex],plotfile,title1,xlabel1,ylabel1,plotmax)
    
plotindex = 0 #[err, stdeverr, L1,L0]
plotmax = np.amax(runArray[:,:,:,plotindex])
x = np.log2(varsList[1])
y = varsList[0] # choose 0,1,2
for i,growvar in enumerate(growlist):
    plotfile = 'growvar%s' % int(100*growvar)
    title1 = 'Validation error for N(o)=%s*N(o-1) , AgPt ' % growvar
    ylabel1 = 'Number of structures in training set'
    xlabel1 = 'Log2 of N-pairs'
    plotArray(x,y,runArray[:,:,i,plotindex],plotfile,title1,xlabel1,ylabel1,plotmax)

################### Sparsity ############
##[structureslist,clusterlist,growlist]
#for ins,ns in enumerate(structureslist):
#    for inc,nc in enumerate(clusterlist):
#        for igw,gw in enumerate(growlist):
#            ntot = nc*(1+gw+gw**2+gw**3+gw**4)
#            spars = ntot/ns
            


print 'Done'





