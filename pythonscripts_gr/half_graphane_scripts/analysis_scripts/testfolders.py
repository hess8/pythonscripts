#print the name of files to analyze
#Specify Directory to use
mainDir = "/bluehome/bch/vasprun/graphene.structures/half_graphane/"

#Specify the subdir
subdir = 'test/'
dir = mainDir + subdir
#Specify the name of the type of run
runName = "relaxation"

#Specify Poscar variables
poscarVariables = {
'@distance':
	[8,6,4,3,2,1]
}

import os,subprocess,time,sys, shutil
 
run = runName
newRun = "DOS"
newRunFile = "DOSCAR" #Will check to see if higher level is already done, then skip it

global toCheckList
global checkedList
global toRunList

def addToList(folder):
    files = os.listdir(folder)
    print files
    for path in files:
#        print path
        if os.path.isdir(folder+path+"/"):
            toCheckList.append(folder+path+"/")
            addToList(folder+path+"/")
#            print path+"/"

def checkFolders():
    for path in toCheckList:
        print("CHECK NEXT LINE")
        print(path.split("/")[-2])
        if path.split("/")[-2] == run:
            checkedList.append(path)

print("\nInitializing...\n")
print("Finding Directories to do %s on\n" % run)
print("Searching all Directories in " + dir+"\n")
toCheckList = []
checkedList = []
toRunList = []
addToList(dir)

print toCheckList
print("Checking to see which folders contain "+run+"\n")
time.sleep(1)
checkFolders()

print "\nThe following folders are in checkedList:"
for i in checkedList:
    print("checkedList contains : " + i)
    
#print "\nThe following folders are in checkedFolder:"
#for i in toCheckList:
#    print("toCheckList contains : " + i)
#print "\nThe following folders will be run:"
#for i in toRunList:
#    print("toRunList contains : " + i)
os.chdir(dir)
file = open("output.txt",'w')
for i in checkedList:
    os.chdir(i)
    print("Testing OSZICAR in: " + i)
    oszicar = open(i+"OSZICAR",'r')
    file.write("\t\t" + oszicar.readlines()[-1].split()[2] + "\n")
    oszicar.close()
file.close()

#Find distance of adatom for minimum energy
os.chdir(dir)
resultsfile = open("output.txt",'r')
results = oszicar.readlines()
for i in range(length

print "Done"
