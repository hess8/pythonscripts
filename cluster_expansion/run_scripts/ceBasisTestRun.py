''' For isolated atoms '''
################## Directories ################## 
#Specify Directories to use
mainDir = "/fslhome/bch/cluster_expansion/cluster_size_test/agpt/"
# input files directory
a = mainDir.split('/')
del a[-2]  #removes last part of path 
inputDir = '/'.join(a) + 'input/'


#Specify the type of run
runtype = 'ncl_ntr'
#runtype = ['test']

#Specify the name of the type of run
runname = 'run_10_15' #make clusters

#Specify a lat.in file
latInFile = inputDir + 'lat.in'

#Specify a CS.in file
csInFile = inputDir + 'CS.in'

################## Variables ################## 

#csInVariables = {
#                  '@NFITSTRUC':[4,8,16,32,64,128,256,512,1024],
#                  '@NFITS':[10]
#                  }

csInVariables = [
                 ['@NFITSTRUC','CS.in',[64,128]],
                 ['@NFITS','CS.in',[10]]
                  ]

#latInVariables = {
#                  '@N2BODY':['CS.in'4,16,64,256,1024],
#                  '@N3BODY':[4,16,64,256,1024],
#                  '@N4BODY':[4,16,64,256,1024],
#                  '@N5BODY':[4,16,64,256,1024],
#                  '@N6BODY':[4,16,64,256,1024],
#                

clusterlist = [4,16]
latInVariables = [
                  ['@N2BODY','lat.in',clusterlist],
                  ['@N3BODY','lat.in',clusterlist],
                  ['@N4BODY','lat.in',clusterlist],
                  ['@N5BODY','lat.in',clusterlist],
                  ['@N6BODY','lat.in',clusterlist]                
                  ]

inputlist = [csInVariables,latInVariables]

################## Build run folders ################## 
import ceScriptTools

toCheckList = []
checkedList = []
toRunList = []
tools = ceScriptTools.ceTools(mainDir,inputDir,runname,runtype,inputlist,toCheckList,checkedList,toRunList)
tools.BuildNewRun() #create folders

raw_input("Press enter to submit jobs")


import os,subprocess,time,sys, shutil

#newRun = "DOS"
#newRunFile = "DOSCAR" #Will check to see if higher level is already done, then skip it
#copyFiles = True


print("\nInitializing...\n")
print("Finding Directories to do %s on\n" % runName)
dir = mainDir + runType[0] + "/"
print("Searching all Directories in " + dir+"\n")
if os.path.exists(dir+'energies'):
    os.system('cp '+dir+'energies '+ dir+'oldenergies')

tools.AddToList(dir)

print("Checking to see which folders contain "+runName+"\n")
time.sleep(1)
tools.CheckFolders()
checkedList = sorted(checkedList)
#toRunList=sorted(toRunList)
toRunList=checkedList #run all folders for convergence checks

#copy INCAR back into relax folder

for folder in toRunList:
    subprocess.call(['cp',incar,folder+'INCAR'])
#    subprocess.call(['rm',folder+'converged.txt'])
#    subprocess.call(['ls','-l',folder+'/'+'INCAR']) 
#    subprocess.call(['cat',folder+'/'+'INCAR'])    
#modify mixing
tools.replaceParamIncar(toRunList, 'BMIX', bmix) # may be needed for convergence in isolated atoms
   
#print "\nThe following folders are in checkedList:"
#for i in checkedList:
#    print("checkedList contains : " + i)
    
print "Skipping the following files, isolated atom converged:"
templist = toRunList[:] #create new list, not just a tag
for i, path in enumerate(toRunList):
    if os.path.exists(path+'converged.txt'):
        templist.remove(path)
        print 'Already converged', path
toRunList = templist[:] 
             
print "\nThe following folders will be run:"
for i in toRunList:
    print("toRunList contains : " + i)

print"\n"
for folder in toRunList:
    newFolder = folder
    os.chdir(newFolder)
    file = open(newFolder+"job",'w+')
    prefix = 'adatom_'
    element = newFolder[newFolder.index(prefix)+len(prefix):newFolder.index('/',newFolder.index(prefix))]
    jobData = "#!/bin/bash\n#PBS -l nodes=1:ppn=1,pmem=2gb,walltime=36:00:00\n#PBS -N " + element+"\n#PBS -m bea\n#PBS -M bret_hess@byu.edu\n# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.\nexport OMP_NUM_THREADS=8\nOUTFILE=\"output.txt\"\n# The following line changes to the directory that you submit your job from\ncd \"$PBS_O_WORKDIR\"\nmpiexec /fslhome/bch/hessgroup/vaspfiles/src/vasp.5.2.12/vasp  > \"$OUTFILE\" \n date > finish.txt \n exit 0"
    file.write(jobData)
    file.close()
    subprocess.check_call(['qsub','job']) #waits to get response 
print "Done with submitting jobs"
