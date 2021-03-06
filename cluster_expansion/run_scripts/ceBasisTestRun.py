''' For isolated atoms '''
import sys,os,subprocess,time,sys,shutil, numpy as np, copy
import ceScriptTools
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
################## Run output to check for completions ################## 
#os.system('python26 /bluehome/bch/pythonscripts/cluster_expansion/analysis_scripts/cdBasisTest_output.py')

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

clusterlist = [4,8,6,16,23,32,45,64,90,128,181,256,362,512,724,1024,1448,2048] #n-2body
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

growlist = [0.0,0.1,0.2,0.3,0.4,0.5,0.56,0.76,1.0,1.3,1.8,2.2]
#growlist = [1.0]
growVariables = [
                  ['*grow','lat.in',growlist]  #for finding n3body, etc...                              
                  ] # use anything but @ for calculated value

structureslist=[32,64,128,256,512,800,900]
#structureslist=[64]
csInVariables = [
                 ['@NFITSTRUC','CS.in',structureslist],
                 ['@NFITS','CS.in',[10]]
                  ]

inputlist = [csInVariables,latInVariables,growVariables]
################## Initialize  ################## 
toCheckList = []
checkedList = []
toRunList = []
tools = ceScriptTools.ceTools(mainDir,inputDir,vaspDataDir,runname,runtype,inputlist,toCheckList,checkedList,toRunList)

os.chdir(dir)
np.save('structureslist',structureslist)
np.save('clusterlist',clusterlist)
np.save('growlist',growlist)
ceScriptTools.writeList(dir,structureslist, 'structureslist.dat')
ceScriptTools.writeList(dir,clusterlist, 'clusterlist.dat')
ceScriptTools.writeList(dir,growlist, 'growlist.dat')
################## Build run folders ################## 
print 'Building run'
tools.BuildNewRun() #create folders




print('Finding Directories to do %s on\n' % runname)

print('Searching all Directories in ' + dir+'\n')

tools.AddToList(dir)

print('Checking to see which folders contain '+runname+'\n')
tools.CheckFolders()
 
toRunList = sorted(tools.RemoveBadRuns(toRunList)) #Removes folders from list and deletes if number of structures is greater than number of clusters
print '\nFolders after remove bad runs:'
for i in toRunList:
    print('toRunList contains : ' + i)
print'\n'
jobname = runtype
outfile = ['out1.txt', 'out2.txt']
execPath = '/fslhome/bch/cluster_expansion/theunclebeta/uncle/trunk/uncle.x'
execLines = '' #initializing
for i,par in enumerate(uncleParam):
    execLines += execPath+' '+str(par)+ ' > ' + outfile[i] + '\n'
print execLines
from ceAnalysisTools import getValues
from ceScriptTools import runningJobs
runninglist = runningJobs()
templist = toRunList[:]
for folder in toRunList:
#    print 'in templist loop', i, folder[-50:]
    nameList = getValues(folder)
    name = ''
    for item in nameList:
        name += str(item)+'_'
#    print name
    if name not in runninglist:        
        os.chdir(folder)    
        file = open(folder+'job','w+')    
        jobData = '#!/bin/bash\n#PBS -l nodes=1:beta,ppn=1,pmem=16gb,walltime=0:20:00\n#PBS -N ' + name +'\n#PBS -m bea\n#PBS -M bret_hess@byu.edu\n' + execLines + 'exit 0'
        file.write(jobData)
        file.close()
    else:
        print 'already running', name
        templist.remove(folder)
toRunList = templist[:]

print '\nThe following folders will be run:'
for i in toRunList:
    print('toRunList contains : ' + i[-50:])
print'\n'

print '%s jobs will be submitted' % len(toRunList)
raw_input('Press enter to submit jobs')
################## Submit jobs ################## 
## remove all ".out files before running"

for folder in toRunList:
    os.chdir(folder)
    print getValues(folder)
    os.system('rm *.out') #so it will write to them later
    #os.system('rm fort.*') 
    subprocess.check_call(['qsub','job'],stderr=subprocess.STDOUT) #waits to get response 
    time.sleep(0.1)
print 'Done with submitting jobs'
