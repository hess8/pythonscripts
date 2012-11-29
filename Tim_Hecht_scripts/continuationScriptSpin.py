#script.py
import os,subprocess,time,sys, shutil
mainDir = "/bluehome/thecht/TransitionMetals/"
previousRun = "IonicRelaxation"
newRun = "DOS"
newRunFile = "DOSCAR"
kpoints = "Automatic mesh\n0\nGamma\n33 33 1\n0 0 0"
incar = "ENCUT=500\nPREC=Accurate\nEDIFF=1E-6\nISMEAR=-5\nSIGMA=0.12\nISPIN=2\n"
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
        if path.split("/")[-2] == previousRun:
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
print("Finding Directories to do DOS run on\n")
print("Searching all Directories in " + mainDir+"\n")
toCheckList = []
checkedList = []
toRunList = []
addToList(mainDir)

print("Checking to see which folders contain "+previousRun+"\n")
time.sleep(1)
checkFolders()
print("Checking to see if any folders have already been run with "+newRun+"\n")
time.sleep(1)
checkForNewRun()
print "\nThe following folders will be run:"
for i in toRunList:
    print i
print"\n"
for folder in toRunList:
    newFolder = folder+newRun+"/"
    print "Copying " + previousRun + " to " + newRun + " in folder " + folder
    command = subprocess.check_call(['cp','-r',folder+previousRun,newFolder])
    print "Copying CONTCAR to POSCAR"
    command = subprocess.check_call(['cp',newFolder+"CONTCAR",newFolder+"POSCAR"])
    print "Setting up KPOINTS file"
    file=open(newFolder+"KPOINTS",'w')
    file.write(kpoints)
    file.close()
    print "Setting up INCAR file"
    file=open(newFolder+"INCAR",'w')
    file.write("System = "+newRun+"\n"+incar)
    file.close()
    file = open(newFolder+"output.txt",'w')
    os.chdir(newFolder)
    proc = subprocess.Popen(['mpirun','-np','8',"vasp.mpi",'>','output.txt'],stdout=subprocess.PIPE)
    while proc.poll() is None:
        output = proc.stdout.readline()
        file.write(output)
        print output,
    file.close()
    

