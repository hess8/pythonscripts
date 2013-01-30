'''After initial relaxation tests, run accurate relaxation at distance of lowest energy. '''

################## Directories ################## 
#Specify Directory to use
mainDir = "/bluehome/bch/vasprun/graphene.structures/transmet.half_graphane/h.half_graphane2.1/"
#mainDir = "/bluehome/bch/vasprun/graphene.structures/transmet.half_graphane/half_graphane/"
#mainDir = "/bluehome/bch/vasprun/graphene.structures/half_graphane/"
#mainDir = "/bluehome/bch/vasprun/graphene.structures/ds_diam_like/"
#vasp input directory
a = mainDir.split('/')

del a[-2]  #removes last part of path 
inputDir = '/'.join(a) + 'vasp.input/'
 
#Specify Potcar Directory
potcardir = "/bluehome/bch/hessgroup/vaspfiles/src/potpaw_PBE/"
#Specify final_relax Directory 
contcardir = mainDir +"final_relax/"

#Specify the type of run
#runType = ['test2']
runType = ['dos']
#Specify the name of run
runName ='dos'

import os,subprocess,time,sys,shutil
if not os.path.exists(mainDir + 'dos/'):
    subprocess.call(['mkdir',mainDir + 'dos/'])
    
#Specify a Poscar file
poscar = inputDir + 'poscar/initialRelax.poscar' #will be replaced by contcar from final relax run

#Specify a KPoints file
kpoints = inputDir + 'kpoints/dos.kpoints'

#Specify an Incar file
incar = inputDir + 'incar/dos.incar'

#Specify a Potcar file
potcar = inputDir + 'potcar/CH.potcar'

################## Variables ################## 

#Specify Poscar variables
poscarVariables = {
#'@distance':
#	[0]
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
'@adatom':
[
'Co','Cr','Cr_pv','Fe','Hf','Ir','Mn','Mn_pv','Mo','Nb_pv','Ni','Os',
'Pd','Pt','Re','Rh','Ru','Ta','Tc','Ti','V','W','Zr'
]
}

#elementList = {
#'@adatom':
#[
#"Os"
#]
#}

################## Build run folders ################## 

import scriptTools
toCheckList = []
checkedList = []
toRunList = []


tools = scriptTools.VaspTools(mainDir,runName,runType,poscar,kpoints,incar,
	potcar,poscarVariables,kpointVariables,incarVariables,elementList,potcardir,
	toCheckList, checkedList, toRunList )

print("\nInitializing...\n")
print("Finding Directories to do %s on\n" % runName)
dir = mainDir + runType[0] + "/"
print("Searching all Directories in " + dir+"\n")
tools.AddToList(dir)

tools.BuildNewRun() #create folders


raw_input("Done creating folders.  Press enter to submit jobs")



################## submit jobs ################## 

#script.py

sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/analysis_scripts')
#mainDir = "/bluehome/bch/TransitionMetals/"
run = runName
newRun = "DOS"
newRunFile = "DOSCAR" #Will check to see if higher level is already done, then skip it
#kpoints = "Automatic mesh\n0\nGamma\n33 33 1\n0 0 0"
#incar = "ENCUT=500\nPREC=Accurate\nEDIFF=1E-6\nISMEAR=-5\nSIGMA=0.12\nISPIN=2\nLORBIT=10\n"
copyFiles = True

print("\nInitializing...\n")
print("Finding Directories to do %s on\n" % run)
print("Searching all Directories in " + dir+"\n")

print("Checking to see which folders contain "+run+"\n")
time.sleep(1)
os.chdir(dir)
print("Searching all Directories in " + dir+"\n")

tools.CheckFolders()
checkedList=sorted(checkedList)
toRunList=sorted(toRunList)

#Copy contcars to poscars
from analysisTools import getElement
pathparts2 = contcardir.split('/')
nparts2 = len(pathparts2)
print nparts2
for path1 in checkedList:
#	element = getElement('adatom_',path)
    pathparts1 = path1.split('/')
#    print pathparts1
    nparts1 = len(pathparts1)#getting num of parts after main dir end
#    print nparts1    
    path2 = contcardir
#    print (pathparts1[nparts2-1:nparts1-1])
    for i, parti in enumerate(pathparts1[nparts2-1:nparts1-2]):
#        print i
#        print path2, parti
        path2 += parti+'/' 
    path2 += 'relax/'
    print("CONTCAR to POSCAR " + path2)
    os.system('cp ' + path2 + 'CONTCAR ' + path1 +'POSCAR')
os.chdir(dir)

print "\nThe following folders will be run:"
for i in toRunList:
    print("toRunList contains : " + i)

print"\n"
for folder in toRunList:
    newFolder = folder
    print folder
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
    prefix = 'adatom_'
    element = newFolder[newFolder.index(prefix)+len(prefix):newFolder.index('/',newFolder.index(prefix))]
    jobData = "#!/bin/bash\n#PBS -l nodes=1:ppn=1,pmem=1gb,walltime=36:00:00\n#PBS -N dos2.1" + element+ "\n#PBS -m bea\n#PBS -M bret_hess@byu.edu\n# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.\nexport OMP_NUM_THREADS=8\nOUTFILE=\"output.txt\"\n# The following line changes to the directory that you submit your job from\ncd \"$PBS_O_WORKDIR\"\nmpiexec /fslhome/bch/hessgroup/vaspfiles/src/vasp.5.2.12/vasp  > \"$OUTFILE\" \n exit 0"
    file.write(jobData)
    file.close()
    subprocess.call(['qsub','job']) #waits to get response 
    

print "Done with submitting jobs"