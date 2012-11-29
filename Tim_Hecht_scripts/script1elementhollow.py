#script.py
import os,subprocess
print("Initializing...\n")
#elementsList = [
#"Ac", "Ca_sv", "Eu", "H.75", "Lu_3", "Np", "Pr_3", "Si_h","Tm",
#"Ac_s", "Cd", "Eu_2", "He", "Mg", "Np_s", "Pt", "Sm", "Tm_3",
#"Ag", "Ce", "F", "Hf", "Mg_pv", "N_s", "Pu", "Sm_3", "U",
#"Al", "Ce_3", "Fe", "Hf_pv", "Mn", "O", "Pu_s", "Sn", "U_s",
#"Al_h", "C_h", "Fe_pv", "Hg", "Mn_pv", "O_h", "Rb_pv", "Sn_d", "V",
#"Ar", "Cl", "F_h", "H_h", "Mo", "Os", "Rb_sv", "Sr_sv", "V_pv",
#"As", "Cl_h", "F_s", "Ho_3", "Mo_pv", "O_s", "Re", "Ta", "V_sv",
#"Au", "Co", "Ga", "I", "N", "Os_pv", "Ta_pv", "W",
#"B", "Cr", "Ga_d", "In", "Na", "P", "Re_pv", "Tb_3", "W_pv",
#"Ba_sv", "Cr_pv", "Ga_h", "In_d", "Na_pv", "Pa", "Rh", "Tc", "Xe",
#"Be", "C_s", "Gd", "Ir", "Na_sv", "Pa_s", "Rh_pv", "Tc_pv", "Yb",
#"Be_sv", "Cs_sv", "Gd_3", "K_pv", "Nb_pv", "Pb", "Ru", "Te", "Yb_2",
#"B_h", "Cu", "Ge", "Kr", "Nb_sv", "Pb_d", "Ru_pv", "Th", "Y_sv",
#"Bi", "Cu_pv", "Ge_d", "K_sv", "Nd", "Pd", "S", "Th_s", "Zn",
#"Bi_d", "Ge_h", "La", "Nd_3", "Pd_pv", "Sb", "Ti", "Zr",
#"Br", "H", "La_s", "Ne", "P_h", "Sc_sv", "Ti_pv", "Zr_sv",
#"B_s", "Dy_3", "H1.25", "Li", "N_h", "Pm", "Se", "Ti_sv",
#"C", "Er_2", "H1.5", "Li_sv", "Ni", "Pm_3", "S_h", "Tl",
#"Ca_pv", "Er_3", "H.5", "Lu", "Ni_pv", "Pr", "Si", "Tl_d"
#]
elementsList = {"B":.87, "N":.56,"H":.53,"C":.67,"O":.48,"Si":.88,
                "F":.42,"Al":1.18,"P":.98,"S":.88,"Cl":.79,"Ge":1.25,
                "Ga":1.36,"As":1.14,"Se":1.03,"Br":.94,"In":1.56,
                "Sn":1.45,"Sb":1.33,"Te":1.23,"I":1.15,"Tl":1.56,"Pb":1.54,
                "Bi":1.43,"Zn":1.42,"Cd":1.61,"Hg":1.71,"Cu":1.45,"Ag":1.65,
                "Au":1.74,"Ni":1.49,"Co":1.52,"Pd":1.69,"Rh":1.73,"Pt":1.77,
                "Ir":1.80,"Fe":1.56,"Ru":1.78,"Os":1.85,"Mn":1.61,"Tc":1.83,
                "Re":1.88,"Cr":1.66,"Mo":1.90,"W":1.93,"V":1.71,
                "Ta":2.00,"Ti":1.76,"Zr":2.06,"Hf":2.08,"Sc_sv":1.84,"Y_sv":2.12,
                "Li":1.67,"Be":1.12,"Na":1.90,"Mg":1.45,"K_sv":2.43,"Ca_sv":1.94,
                "Sr_sv":2.19,"Cs_sv":2.98,"Rb_sv":2.35,"Ba_sv":2.15,"La":2.12,
                "Nb_pv":1.98}
                
outputFile = open("output1elementhollow.dat",'w')

totalFiles=len(elementsList)*2
print str(len(elementsList)*2) + " files to be run"
curFile=1
poscar1 = "1\n2.13129 -1.2305 0\n2.13129 1.2305 0\n0.0 0.0 30.0\n"
poscar2 = "Cartesian\n1.42086 0 0\n2.84172 0 0\n"
incar = "IBRION=2\nISIF=4\nPOTIM=0.5\nNSW=50\nENCUT=500\nPREC=Accurate\nEDIFF=1E-6\nISMEAR=1"
kpoints = "Automatic mesh\n0\nGamma\n11 11 1\n0 0 0"
try:
    os.chdir("1elementhollow")
except OSError as e:
    os.mkdir("1elementhollow")
    os.chdir("1elementhollow")
    
def checkFolder(folder):
    global curFile
    print("Starting " + folder +"...")
    print("File " + str(curFile) + " out of " + str(totalFiles))
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

def Poscar(isa1,isa2):
    file = open("POSCAR", 'w')
    text = folder + "\n"+poscar1
    num = 0
    temp = poscar2
    dist = .7 + elementsList[label]
    if isa1:
        num = num+1
        temp=temp+"0 0 "+str(dist)+"\n"
    if isa2:
        num= num+1
        temp=temp+"0 0 "+str(-dist)+"\n"
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
    os.chdir("/bluehome/thecht/highoutput/1elementhollow")
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
    file = open("output.txt",'w')
    proc = subprocess.Popen(['mpirun','-np','8','vasp.mpi','>','output.txt'],stdout=subprocess.PIPE)
    while proc.poll() is None:
        output = proc.stdout.readline()
        file.write(output)
        print output,
    oszicar = open("OSZICAR",'r')
    outputFile.write(folder + "\t\t" + oszicar.readlines()[-1].split()[2] + "\n")
    oszicar.close()
    file.close()


def setup(folder,b1,b2):
    if checkFolder(folder)==True:
        print("Preparing Files...")
        Poscar(b1,b2)
        Incar()
        Potcar()
        Kpoints()
        SymbolicLink()
        run()

for label in elementsList.keys():
    folder = label + ":"
    setup(folder,True,False)
    os.chdir("..")
    folder = label + ":" + label
    setup(folder,True,True)
    os.chdir("..")


outputFile.close()
