#script.py
import os,subprocess
print("Initializing...\n")
elementsList = [
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
outputFile = open("outputstretch.dat",'w')

print str(len(elementsList)*6) + " files to be run"
numFiles = len(elementsList)*6
curFile = 1
poscar1 = "1\n2.13129 -1.2305 0\n2.13129 1.2305 0\n0.0 0.0 30.0\n"
poscar2 = "Cartesian\n1.42086 0 0\n2.84172 0 0\n"
incar = "IBRION=-1\nNELM=100\nEDIFF = 1E-06\nPREC=Accurate\nENCUT=500\n"
kpoints = "Automatic mesh\n0\nGamma\n11 11 1\n0 0 0"
try:
    os.chdir("stretch")
except OSError as e:
    os.mkdir("stretch")
    os.chdir("stretch")
    
def checkFolder(folder):
    global curFile
    print("Starting " + folder +"...")
    print(str(curFile) + " out of " + str(numFiles))
    curFile = curFile + 1
    try:
       os.chdir(folder)
       print("Skipping " + folder +", folder was already found.\n")
       oszicar = open("OSZICAR",'r')
       outputFile.write(folder + "\t\t" + oszicar.readlines()[-1].split()[2] + "\n")
       oszicar.close()
       return False
    except OSError as e:
        os.mkdir(folder)
        os.chdir(folder)
        return True

def Poscar(isa1,isb1,isa2,isb2):
    file = open("POSCAR", 'w')
    text = folder + "\n"+poscar1
    num = 0
    temp = poscar2
    dist = 8.00
    if isa1:
        num = num+1
        temp=temp+"1.42086 0 "+str(dist)+"\n"
    if isb1:
        num= num+1
        temp=temp+"2.84172 0 "+str(dist)+"\n"
    if isa2:
        num= num+1
        temp=temp+"1.42086 0 "+str(-dist)+"\n"
    if isb2:
        num= num+1
        temp=temp+"2.84172 0 "+str(-dist)+"\n"
    text2 = text + "2 " + str(num) + "\n" + temp
    file.write(text2)
    file.close()

def Incar():
    file=open("INCAR",'w')
    text = "System = "+folder+"\n"+incar
    file.write(text)
    file.close()

def Potcar():
    os.chdir("/bluehome/thecht/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/")
    os.chdir(label)
    file=open("POTCAR",'r')
    file2=file.readlines()
    file.close()
    file2string=""
    for string in file2:
        file2string=file2string+string
    os.chdir("..")
    os.chdir("C")
    file=open("POTCAR",'r')
    file1=file.readlines()
    file.close()
    file1string=""
    for string in file1:
        file1string=file1string+string
    os.chdir("/bluehome/thecht/highoutput/stretch")
    os.chdir(folder)
    file=open("POTCAR",'w')
    file.write(file1string+file2string)
    file.close()

def Kpoints():
    file=open("KPOINTS",'w')
    text = kpoints
    file.write(text)
    file.close()


def SymbolicLink():
    p = subprocess.check_call(['cp','.././../vasp.mpi','.'])


def run():
    proc = subprocess.Popen(['mpirun','-np','8','vasp.mpi'],stdout=subprocess.PIPE)
    while proc.poll() is None:
        output = proc.stdout.readline()
        print output,
    oszicar = open("OSZICAR",'r')
    outputFile.write(folder + "\t\t" + oszicar.readlines()[-1].split()[2] + "\n")
    oszicar.close()


def setup(folder,b1,b2,b3,b4):
    if checkFolder(folder)==True:
        print("Preparing Files...")
        Poscar(b1,b2,b3,b4)
        Incar()
        Potcar()
        Kpoints()
        SymbolicLink()
        run()

for label in elementsList:
    folder = label + ":"
    setup(folder,True,False,False,False)
    os.chdir("..")
    folder = label + "-" + label + ":"
    setup(folder,True,True,False,False)
    os.chdir("..")
    folder = label + ":" + label
    setup(folder,True,False,True,False)
    os.chdir("..")
    folder = label + ":-" + label
    setup(folder,True,False,False,True)
    os.chdir("..")
    folder = label + "-" + label + ":" + label 
    setup(folder,True,True,True,False)
    os.chdir("..")
    folder = label + "-" + label + ":" + label + "-" + label
    setup(folder,True,True,True,True)
    os.chdir("..")

outputFile.close()

