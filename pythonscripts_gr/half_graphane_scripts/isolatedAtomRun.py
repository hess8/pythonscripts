''' For isolated atoms '''
################## Directories ################## 
#Specify Directory to use
mainDir = "/bluehome/bch/vasprun/graphene.structures/transmet.half_graphane/isolated/"
a = mainDir.split('/')
del a[-2]  #removes last part of path 
inputDir = '/'.join(a) 
print inputDir
bmix = '0.001'
#Specify Potcar Directory
potcardir = "/bluehome/bch/hessgroup/vaspfiles/src/potpaw_PBE/"

#Specify the type of run
runType = ['electron_relax']
#runType = ['test']

#Specify the name of the type of run
runName = "relaxation"

#Specify a Poscar file
poscar = inputDir + 'vasp.input/poscar/isolatedrun.poscar'

#Specify a KPoints file
kpoints = inputDir + 'vasp.input/kpoints/isolatedrun.kpoints'

#Specify an Incar file
incar = inputDir + 'vasp.input/incar/isolatedrun.incar'

#Specify a Potcar file
potcar = inputDir + 'vasp.input/potcar/isolatedrun.potcar'

################## Variables ################## 

#Specify Poscar variables
poscarVariables = {
}

#Specify KPoints variables
kpointVariables = {
}

#Specify Incar variables
incarVariables = {
}



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
	toCheckList, checkedList, toRunList)
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
