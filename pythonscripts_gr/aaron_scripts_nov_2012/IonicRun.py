#Specify Directory to use
mainDir = "/bluehome/aaronm3/TransitionMetals/"

#Specify Potcar Directory
potcardir = "/bluehome/aaronm3/hessgroup/vaspfiles/src/potpaw_PBE/"

#Specify the type of run
runType = ['Relax']

#Specify the name of the type of run
runName = "Relaxation"

#Specify a Poscar file
poscar = mainDir + 'Poscar/Relax.poscar'

#Specify a KPoints file
kpoints = mainDir + 'Kpoints/grapheneionicrun.kpoints'

#Specify an Incar file
incar = mainDir + 'Incar/grapheneionicrun.incar'

#Specify a Potcar file
potcar = mainDir + 'Potcar/graphene.potcar'

#Specify Poscar variables
poscarVariables = {
'@distance':
	[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0]
}

#Specify KPoints variables
kpointVariables = {

}

#Specify Incar variables
incarVariables = {
'@ISPIN':
	['2'],
'@IBRION':
	['2'],
'@ISIF':
	['4']

}

#Specify Potcar Elements
elementList = {
'@Adatom':
	['Hf','La']
}

import ScriptTools

tools = ScriptTools.VaspTools(mainDir,runName,runType,poscar,kpoints,incar,
	potcar,poscarVariables,kpointVariables,incarVariables,elementList,potcardir)

tools.BuildNewRun()

print "Done"
raw_input()












#script.py
import os,subprocess,time,sys, shutil
mainDir = "/bluehome/aaronm3/TransitionMetals/"
run = runName
newRun = "DOS"
newRunFile = "DOSCAR"
kpoints = "Automatic mesh\n0\nGamma\n33 33 1\n0 0 0"
incar = "ENCUT=500\nPREC=Accurate\nEDIFF=1E-6\nISMEAR=-5\nSIGMA=0.12\nISPIN=2\nLORBIT=10\n"
copyFiles = True

global toCheckList
global checkedList
global toRunList

def addToList(folder):
    files = os.listdir(folder)
    for path in files:
        if os.path.isdir(folder+path+"/"):
            toCheckList.append(folder+path+"/")
            addToList(folder+path+"/")

def checkFolders():
    for path in toCheckList:
        print("CHECK NEXT LINE")
        print(path.split("/")[-2])
        if path.split("/")[-2] == run:
            checkedList.append(path)
            
def checkForNewRun():
    for path in checkedList:
        parpath =  os.path.abspath(os.path.join(path, os.path.pardir))
        if os.path.exists(os.path.join(parpath,newRun)):
            print os.path.join(parpath,newRun) + " already exists."
            if copyFiles:
                print "Copying " + newRunFile +" from path to current directory."
                newPath = os.path.join(parpath,newRun) + "/" + newRunFile
                array = parpath.split("/")
                newFileName = array[len(array)-2]+array[len(array)-1]+".dat"
                shutil.copyfile(newPath,mainDir+newFileName)
        else:
            toRunList.append(parpath+"/")

print("\nInitializing...\n")
print("Finding Directories to do %s on\n" % run)
print("Searching all Directories in " + mainDir+"\n")
toCheckList = []
checkedList = []
toRunList = []
addToList(mainDir)

print("Checking to see which folders contain "+run+"\n")
time.sleep(1)
checkFolders()

print("Checking to see if any folders have already converged\n")
time.sleep(1)
checkForNewRun()

print "\nThe following folders are in checkedList:"
for i in checkedList:
    print("checkedList contains : " + i)
    
print "\nThe following folders are in checkedFolder:"
for i in toCheckList:
    print("toCheckList contains : " + i)
print "\nThe following folders will be run:"
for i in toRunList:
    print("toRunList contains : " + i)

print "\nThe script is at line 134\n"

print"\n"
for folder in toRunList:
    newFolder = folder+run+"/"
    #print "Copying " + previousRun + " to " + newRun + " in folder " + folder
    #command = subprocess.check_call(['cp','-r',folder+previousRun,newFolder])
    #print "Copying CONTCAR to POSCAR"
    #command = subprocess.check_call(['cp',newFolder+"CONTCAR",newFolder+"POSCAR"])
    #print "Setting up KPOINTS file"
    #file=open(newFolder+"KPOINTS",'w')
    #file.write(kpoints)
    #file.close()
    #print "Setting up INCAR file"
    #file=open(newFolder+"INCAR",'w')
    #file.write("System = "+newRun+"\n"+incar)
    #file.close()
    os.chdir(newFolder)
    file = open(newFolder+"job",'w+')
    #print "made job in"
    #print newFolder + "job"
    jobData = "#!/bin/bash\n#PBS -l nodes=1:ppn=8,pmem=60mb,walltime=00:5:00\n#PBS -N JOBNAME\n#PBS -m bea\n#PBS -M l3370x@gmail.com\n# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.\nexport OMP_NUM_THREADS=8\nOUTFILE=\"output.txt\"\n# The following line changes to the directory that you submit your job from\ncd \"$PBS_O_WORKDIR\"\nmpiexec /fslhome/aaronm3/vasp.mpi  > \"$OUTFILE\" \nexit 0"
    file.write(jobData)
    file.close()
    file = open(newFolder+"outputJob.txt",'w')
    proc = subprocess.Popen(['qsub','job'],stdout=subprocess.PIPE)
    while proc.poll() is None:
        output = proc.stdout.readline()
        file.write(output)
        print output,folder
    file.close()
    

