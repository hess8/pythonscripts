''' For isolated atoms '''
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
#clusterlist = [4,8,16,32] #,64,128,256,512,1024,2048] 
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

growVariables = [
                  ['*grow','lat.in',[0.56,0.76,1.0,1.3,1.8,]]  #for finding n3body, etc...                              
                  ] # use anything but @ for calculated value


csInVariables = [
                 ['@NFITSTRUC','CS.in',[64,128,256,512]],
                 ['@NFITS','CS.in',[10]]
                  ]

inputlist = [csInVariables,latInVariables,growVariables]

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
dir = mainDir + runtype + '/'
print('Searching all Directories in ' + dir+'\n')

tools.AddToList(dir)

print('Checking to see which folders contain '+runname+'\n')
tools.CheckFolders()
#checkedList = sorted(checkedList)
#toRunList=checkedList  
toRunList = tools.RemoveBadRuns(checkedList) #Removes folders from list and deletes if number of structures is greater than number of clusters
print '\nThe following folders are in checkedList:'
for i in checkedList:
    print('checkedList contains : ' + i)
    
#print 'Skipping the following files, found finish.txt:'
#templist = toRunList[:] #create new list, not just a tag
#for i, path in enumerate(toRunList):
#    if os.path.exists(path+'finish.txt'):
#        templist.remove(path)
#        print 'Already finished', path
#toRunList = templist[:] 
             
print '\nThe following folders will be run:'
for i in toRunList:
    print('toRunList contains : ' + i)
print'\n'
jobname = runtype
outfile = 'out.txt'
execPath = '/fslhome/bch/cluster_expansion/theuncle/uncle/trunk/uncle.x'
execLines = ''
for par in uncleParam:
    execLines += 'mpiexec '+execPath+' '+str(par)+ ' > ' + outfile + '\n'
print execLines
for folder in toRunList:
    os.chdir(folder)
    file = open(folder+'job','w+')    
    jobData = '#!/bin/bash\n#PBS -l nodes=1:ppn=1,pmem=4gb,walltime=4:00:00\n#PBS -N ' + jobname +'\n#PBS -m bea\n#PBS -M bret_hess@byu.edu\n' + execLines + 'date > finish.txt \n exit 0'
    file.write(jobData)
    file.close()
print '%s jobs will be submited' % len(toRunList)
raw_input('Press enter to submit jobs')

for folder in toRunList:
    os.chdir(folder)
    subprocess.check_call(['qsub','job']) #waits to get response 
print 'Done with submitting jobs'
