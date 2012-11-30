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

elementsList = ["B", "N","H","C","O","Si","F","Al","P","S","Cl","Ge",
                "Ga","As","Se","Br","In","Sn","Sb","Te","I","Tl","Pb",
                "Bi","Zn","Cd","Hg","Cu","Ag","Au","Ni","Co","Pd","Rh","Pt",
                "Ir","Fe","Ru","Os","Mn","Tc","Re","Cr","Mo","W","V",
                "Ta","Ti","Zr","Hf","Sc_sv","Y_sv","Li","Be","Na","Mg","K_sv","Ca_sv",
                "Sr_sv","Cs_sv","Rb_sv","Ba_sv","La","Nb_pv"]

outputFile = open("outputexfoliation.dat",'w')

print str(2*len(elementsList)) + " files to be run"
numFiles =2* len(elementsList)
curFile = 1
poscar1 = "1\n2.13129 -1.2305 0\n2.13129 1.2305 0\n0.0 0.0 30.0\n1\nCartesian\n0 0 0\n"
poscar2 = "1\n2.13129 -1.2305 0\n2.13129 1.2305 0\n0.0 0.0 30.0\n2\nCartesian\n1.42086 0 0\n3.2 0 0\n"
incar = "NSW=100\nIBRION=2\nISIF=4\nNELM=100\nEDIFF = 1E-06\nPREC=Accurate\nENCUT=500\n"
kpoints = "Automatic mesh\n0\nGamma\n11 11 1\n0 0 0"
try:
    os.chdir("exfoliation")
except OSError as e:
    os.mkdir("exfoliation")
    os.chdir("exfoliation")
    
def checkFolder(folder):
    global curFile
    print("Starting " + folder +"...")
    print(str(curFile) + " out of " + str(numFiles))
    curFile = curFile + 1
    try:
       os.chdir(folder)
       print("Skipping " + folder +", folder was already found.\n")
       oszicar = open("OSZICAR",'r')
       outputFile.write(folder + "\t\t" + oszicar.readlines()[-1].split()[2]+"\t\t")
       oszicar.close()
       proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
       newstring = proc.communicate()
       numions = int(newstring[0].split()[-1])
       proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
       line = proc2.communicate()[-2].split('\n')[-2].split(':')[0]
       outcar = open('OUTCAR','r')
       text = outcar.readlines()
       outcar.close()
       maximum=float(text[int(line)+1].split()[5])
       outputFile.write(str(maximum)+'\n')
       return False
    except OSError as e:
        os.mkdir(folder)
        os.chdir(folder)
        return True

def Poscar(i):
    file = open("POSCAR", 'w')
    if i ==1:
        text = folder + "\n"+poscar1
    if i ==2:
        text = folder + "\n"+poscar2
    file.write(text)
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
    os.chdir("/bluehome/thecht/highoutput/exfoliation")
    os.chdir(folder)
    file=open("POTCAR",'w')
    file.write(file2string)
    file.close()

def Kpoints():
    file=open("KPOINTS",'w')
    text = kpoints
    file.write(text)
    file.close()


def SymbolicLink():
    p = subprocess.check_call(['ln','-s','/bluehome/thecht/TransitionMetals/vasp.mpi','vasp.mpi'])


def run():
    proc = subprocess.Popen(['mpirun','-np','8','vasp.mpi'],stdout=subprocess.PIPE)
    while proc.poll() is None:
        output = proc.stdout.readline()
        print output,
    oszicar = open("OSZICAR",'r')
    outputFile.write(folder + "\t\t" + oszicar.readlines()[-1].split()[2] + "\t\t")
    proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    numions = int(newstring[0].split()[-1])
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    line = proc2.communicate()[-2].split('\n')[-2].split(':')[0]
    outcar = open('OUTCAR','r')
    text = outcar.readlines()
    outcar.close()
    maximum=float(text[int(line)+1].split()[5])
    for i in range(2,numions+1):
      temp=float(text[int(line)+i].split()[5])
      maximum=max(abs(maximum),abs(temp))
    outputFile.write(str(maximum)+'\n')
    oszicar.close()


def setup(folder,i):
    if checkFolder(folder)==True:
        print("Preparing Files...")
        Poscar(i)
        Incar()
        Potcar()
        Kpoints()
        SymbolicLink()
        run()

for label in elementsList:
    folder = label + ":"
    setup(folder,1)
    os.chdir("..")
    folder = label + "-" + label + ":"
    setup(folder,2)
    os.chdir("..")

outputFile.close()

