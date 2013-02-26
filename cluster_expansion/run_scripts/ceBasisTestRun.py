''' For isolated atoms '''
import sys, os
sys.path.append("/fslhome/bch/pythonscripts/cluster_expansion/analysis_scripts/")
################## Directories ################## 
#Specify Directories to use
mainDir = '/fslhome/bch/cluster_expansion/cluster_size_test/agpt/'
# input files directory
a = mainDir.split('/')
del a[-2]  #removes last part of path 
inputDir = '/'.join(a) + 'input/'
vaspDataDir = mainDir+'vasp_data/'

#Specify the type of run
runtype = 'ncl_ntr'
dir = mainDir + runtype + '/'
#Specify the name of the type of run
uncleParam = [10,15]
runname = 'run'
for par in uncleParam:
    runname += '_%s' % par

#Specify a lat.in file
latInFile = inputDir + 'lat.in'

#Specify a CS.in file
csInFile = inputDir + 'CS.in'

################## Variables ################## 

#                

clusterlist = [4,8,16,32,64,128,256,512,1024,2048] 
#clusterlist = [64,128] 

#latInVariables = [
#                  ['@N2BODY','lat.in',clusterlist],
#                  ['@N3BODY','lat.in',clusterlist],
#                  ['@N4BODY','lat.in',clusterlist],
#                  ['@N5BODY','lat.in',clusterlist],
#                  ['@N6BODY','lat.in',clusterlist]                
#                  ]

latInVariables = [
                  ['@N2BODY','lat.in',clusterlist] #will calculate others from growVars                              
                  ]

growlist = [0.56,0.76,1.0,1.3,1.8]
#growlist = [1.0]
growVariables = [
                  ['*grow','lat.in',growlist]  #for finding n3body, etc...                              
                  ] # use anything but @ for calculated value

structureslist=[64,128,256,512,800]
#structureslist=[64]
csInVariables = [
                 ['@NFITSTRUC','CS.in',structureslist],
                 ['@NFITS','CS.in',[10]]
                  ]

inputlist = [csInVariables,latInVariables,growVariables]

os.chdir(dir)
varsFile = open('variables.out','w')
print [round(i,2) for i in growlist]
varsFile.writelines(str([structureslist,clusterlist,[round(i,2) for i in growlist]]))  #write variables list which will be read during analysis
varsFile.close()
################## Build run folders ################## 
import ceScriptTools

toCheckList = []
checkedList = []
toRunList = []
tools = ceScriptTools.ceTools(mainDir,inputDir,vaspDataDir,runname,runtype,inputlist,toCheckList,checkedList,toRunList)
print 'Building run'
tools.BuildNewRun() #create folders


import os,subprocess,time,sys, shutil

print('Finding Directories to do %s on\n' % runname)

print('Searching all Directories in ' + dir+'\n')

tools.AddToList(dir)

print('Checking to see which folders contain '+runname+'\n')
tools.CheckFolders()
 
toRunList = tools.RemoveBadRuns(toRunList) #Removes folders from list and deletes if number of structures is greater than number of clusters
#print '\nThe following folders are in checkedList:'
#for i in checkedList:
#    print('checkedList contains : ' + i)
    
#print 'Skipping the following files, found finish.txt:'
#templist = toRunList[:] #create new list, not just a tag
#for i, path in enumerate(toRunList):
#    if os.path.exists(path+'finish.txt'):
#        templist.remove(path)
#        print 'Already finished', path
#toRunList = templist[:] 
################## Prepare jobs ##################              
print '\nThe following folders will be run:'
for i in toRunList:
    print('toRunList contains : ' + i)
print'\n'
jobname = runtype
outfile = ['out1.txt', 'out2.txt']
execPath = '/fslhome/bch/cluster_expansion/theunclebeta/uncle/trunk/uncle.x'
#execLines = 'export LD_LIBRARY_PATH=/opt/intel/mkl/10.2.5.035/lib/em64t:$LD_LIBRARY_PATH\n' #initializing
execLines = '' #initializing
for i,par in enumerate(uncleParam):
    execLines += execPath+' '+str(par)+ ' > ' + outfile[i] + '\n'
print execLines
from ceAnalysisTools import getValues
for folder in toRunList:
    nameList = getValues(folder)
    name = ''
    for item in nameList:
        name += str(item)+'_'
    print name
    os.chdir(folder)    
    file = open(folder+'job','w+')    
    jobData = '#!/bin/bash\n#PBS -l nodes=1:beta,ppn=1,pmem=16gb,walltime=36:00:00\n#PBS -N ' + name +'\n#PBS -m bea\n#PBS -M bret_hess@byu.edu\n' + execLines + 'exit 0'
    file.write(jobData)
    file.close()
print '%s jobs will be submitted' % len(toRunList)
raw_input('Press enter to submit jobs')
################## Submit jobs ################## 
for folder in toRunList:
    os.chdir(folder)
    subprocess.check_call(['qsub','job']) #waits to get response 
    time.sleep(0.1)
print 'Done with submitting jobs'
