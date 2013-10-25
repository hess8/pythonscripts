'''After initial relaxation tests, run accurate relaxation at distance of lowest energy. '''

################## Directories ################## 
#Specify Directory to use
mainDir = "/fslhome/bch/vasprun/graphene.structures/transmet.half_graphane/h.half_graphane2.1/"
#mainDir = "/bluehome/bch/vasprun/graphene.structures/transmet.half_graphane/half_graphane/"
#mainDir = "/bluehome/bch/vasprun/graphene.structures/half_graphane/"
#mainDir = "/bluehome/bch/vasprun/graphene.structures/ds_diam_like/"


contcardir = '/fslhome/bch/vasprun/graphene.structures/transmet.half_graphane/h.half_graphane2.1/final_relax/adatom_Co/ISIF_4/IBRION_2/ISPIN_2/relax/'
#contcardir = '/fslhome/bch/vasprun/graphene.structures/transmet.half_graphane/half_graphane/final_relax/adatom_Co/ISIF_4/IBRION_2/ISPIN_2/relax/'

#vasp input directory
a = mainDir.split('/')

del a[-2]  #removes last part of path 
inputDir = '/'.join(a) + 'vasp.input/'
 
#Specify Potcar Directory
potcardir = "/fslhome/bch/hessgroup/vaspfiles/src/potpaw_PBE/"
#Specify sample directory where contcar is for replacement of element

#Specify the type of run
#runType = ['test']
runType = ['layer']
#Specify the name of run
runName ='relax'

import os,subprocess,time,sys,shutil
if not os.path.exists(mainDir + runType[0]):
    subprocess.call(['mkdir',mainDir + runType[0]])
    
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
	['2']
}

#Specify Potcar Elements

elementList = {
'@adatom':
[
'Co','Cr','Fe','Hf','Ir'
,'Mn','Mo','Nb_pv','Ni','Os',
'Pd','Pt','Re','Rh','Ru','Ta','Tc','Ti','V','W','Zr'
]
}



#elementList = {
#'@adatom':
#[
#"Cr"
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


#raw_input("Done creating folders.  Press enter to submit jobs")



################## submit jobs ################## 

#script.py

sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/analysis_scripts')
#mainDir = "/bluehome/bch/TransitionMetals/"
run = runName

print("\nInitializing...\n")
print("Finding Directories to do %s on\n" % run)
print("Searching all Directories in " + dir+"\n")

print("Checking to see which folders contain "+run+"\n")
time.sleep(1)
os.chdir(dir)
print("Searching all Directories in " + dir+"\n")

tools.CheckFoldersDOS()
checkedList=sorted(checkedList)
toRunList=sorted(toRunList)

#Copy contcars to poscars
from analysisTools import getElement
#pathparts2 = contcardir.split('/')
#nparts2 = len(pathparts2)
#print nparts2
for path1 in checkedList:
    path2 = contcardir.replace('Co',getElement('adatom_',path1))
    print("cp CONTCAR" + path2)
    print('to POSCAR'+ path1)
    os.system('cp ' + path2 + 'CONTCAR ' + path1 +'POSCAR')
os.chdir(dir)

tools.removeCsPOSCAR(toRunList)

print "\nThe following folders will be run:"
for i in toRunList:
    print("toRunList contains : " + i)
raw_input("Press enter to submit jobs")
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
    jobData = "#!/bin/bash\n#PBS -l nodes=1:beta,ppn=1,pmem=1gb,walltime=36:00:00\n#PBS -N stretch_h." + element+ "\n#PBS -m bea\n#PBS -M bret_hess@byu.edu\n# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.\nexport OMP_NUM_THREADS=8\nOUTFILE=\"output.txt\"\n# The following line changes to the directory that you submit your job from\ncd \"$PBS_O_WORKDIR\"\nmpiexec /fslhome/bch/hessgroup/vaspfiles/src/vasp.5.3.3/vasp  > \"$OUTFILE\" \n exit 0"
    file.write(jobData)
    file.close()
    subprocess.call(['qsub','job']) #waits to get response 
    

print "Done submitting jobs"