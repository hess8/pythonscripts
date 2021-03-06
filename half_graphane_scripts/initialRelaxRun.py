#Specify Directory to use
mainDir = "/bluehome/bch/vasprun/graphene.structures/ds_diam_like/"
#Specify the subdir
runType  = ['initial_relax']

#Specify Potcar Directory
potcardir = "/bluehome/bch/hessgroup/vaspfiles/src/potpaw_PBE/"

#Specify the name of the type of run
runName = 'relax' 

#Specify a Poscar file
poscar = mainDir + 'poscar/initialRelax.poscar'

#Specify a KPoints file
kpoints = mainDir + 'kpoints/initialRelax.kpoints'

#Specify an Incar file
incar = mainDir + 'incar/initialRelax.incar'

#Specify a Potcar file
potcar = mainDir + 'potcar/C.potcar'

#Specify Poscar variables
poscarVariables = {
'@distance':
	[8,6,4,3,2,1]
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

elementList = {
'@adatom':
[
"Ac", "Ca_sv", "Eu", "H.75", "Lu_3", "Np", "Pr_3", "Si_h","Tm",
"Ac_s", "Cd", "Eu_2", "He", "Mg", "Np_s", "Pt", "Sm", "Tm_3",
"Ag", "Ce", "F", "Hf", "Mg_pv", "N_s", "Pu", "Sm_3", "U",
"Al", "Ce_3", "Fe", "Hf_pv", "Mn", "O", "Pu_s", "Sn", "U_s",
"Al_h", "C_h", "Fe_pv", "Hg", "Mn_pv", "O_h", "Rb_pv", "Sn_d", "V",
"Ar", "Cl", "F_h", "H_h", "Mo", "Os", "Rb_sv", "Sr_sv", "V_pv",
"As", "Cl_h", "F_s", "Ho_3", "Mo_pv", "O_s", "Re", "Ta", "V_sv",
"Au", "Co", "Ga", "I", "N", "Os_pv", "Ta_pv", "W",
"B", "Cr", "Ga_d", "In", "Na", "P", "Re_pv", "Tb_3", "W_pv",
"Ba_sv", "Cr_pv", "Ga_h", "In_d", "Na_pv", "Pa", "Rh", "Tc", "Xe",
"Be", "C_s", "Gd", "Ir", "Na_sv", "Pa_s", "Rh_pv", "Tc_pv", "Yb",
"Be_sv", "Cs_sv", "Gd_3", "K_pv", "Nb_pv", "Pb", "Ru", "Te", "Yb_2",
"B_h", "Cu", "Ge", "Kr", "Nb_sv", "Pb_d", "Ru_pv", "Th", "Y_sv",
"Bi", "Cu_pv", "Ge_d", "K_sv", "Nd", "Pd", "S", "Th_s", "Zn",
"Bi_d", "Ge_h", "La", "Nd_3", "Pd_pv", "Sb", "Ti", "Zr",
"Br", "H", "La_s", "Ne", "P_h", "Sc_sv", "Ti_pv", "Zr_sv",
"B_s", "Dy_3", "H1.25", "Li", "N_h", "Pm", "Se", "Ti_sv",
"C", "Er_2", "H1.5", "Li_sv", "Ni", "Pm_3", "S_h", "Tl",
"Ca_pv", "Er_3", "H.5", "Lu", "Ni_pv", "Pr", "Si", "Tl_d"
]
}

#elementList = {
#'@adatom':
#[
#"Ac", "Ca_sv"
#]
#}

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

#script.py
import os,subprocess,time,sys, shutil
#mainDir = "/bluehome/bch/TransitionMetals/"
run = runName
newRun = "DOS"
newRunFile = "DOSCAR" #Will check to see if higher level is already done, then skip it
#kpoints = "Automatic mesh\n0\nGamma\n33 33 1\n0 0 0"
#incar = "ENCUT=500\nPREC=Accurate\nEDIFF=1E-6\nISMEAR=-5\nSIGMA=0.12\nISPIN=2\nLORBIT=10\n"
copyFiles = True


print("\nInitializing...\n")
print("Finding Directories to do %s on\n" % run)
print("Searching all Directories in " + mainDir+"\n")
tools.AddToList(mainDir)

print("Checking to see which folders contain "+run+"\n")
time.sleep(1)
tools.CheckFolders()
checkedList=sorted(checkedList)
toRunList=sorted(toRunList)
print("Checking to see if any folders have already converged\n")
#time.sleep(1)
#checkForNewRun()

print "\nThe following folders are in checkedList:"
for i in checkedList:
    print("checkedList contains : " + i)
    
print "\nThe following folders are in checkedFolder:"
for i in toCheckList:
    print("toCheckList contains : " + i)
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
    #print "made job in"
    #print newFolder + "job"
    prefix = 'adatom_'
    element = newFolder[newFolder.index(prefix)+len(prefix):newFolder.index('/',newFolder.index(prefix))]
    print element
    jobData = "#!/bin/bash\n#PBS -l nodes=1:beta,ppn=1,pmem=1gb,walltime=06:00:00\n#PBS -N fin_rel_" + element+ "\n#PBS -m bea\n#PBS -M bret_hess@byu.edu\n# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.\nexport OMP_NUM_THREADS=8\nOUTFILE=\"output.txt\"\n# The following line changes to the directory that you submit your job from\ncd \"$PBS_O_WORKDIR\"\nmpiexec /fslhome/bch/hessgroup/vaspfiles/src/vasp.5.3.3/vasp  > \"$OUTFILE\" \nexit 0"
    file.write(jobData)
    file.close()
    file = open(newFolder+"outputJob.txt",'w')
    proc = subprocess.Popen(['qsub','job'],stdout=subprocess.PIPE)
    while proc.poll() is None:
        output = proc.stdout.readline()
        file.write(output)
        print output,folder
    file.close()
    
print "Done submitting jobs"
