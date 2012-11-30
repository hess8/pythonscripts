#script.py
import os,subprocess,math
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
file = open('outputcarbon.dat','w')
maindir = "/fslhome/thecht/highoutput/1elementhollow"
os.chdir(maindir)

def run():
    file.write(folder)
    print folder
    proc = subprocess.Popen(['grep','-i','nion','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc.communicate()
    numions = int(newstring[0].split()[-1])
    proc3 = subprocess.Popen(['grep','-i','A3','OUTCAR'],stdout=subprocess.PIPE)
    newstring = proc3.communicate()
    repeat = math.fabs(float(newstring[0].split()[-1].split(")")[0]))
    proc2 = subprocess.Popen(['grep','-n','TOTAL-FORCE (eV/Angst)','OUTCAR'],stdout=subprocess.PIPE)
    line = proc2.communicate()[-2].split('\n')[-2].split(':')[0]
    outcar = open('OUTCAR','r')
    text = outcar.readlines()
    outcar.close()
    carbon1=[float(text[int(line)+1].split()[0]),float(text[int(line)+1].split()[1]),float(text[int(line)+1].split()[2])]
    carbon2=[float(text[int(line)+2].split()[0]),float(text[int(line)+2].split()[1]),float(text[int(line)+2].split()[2])]
    distance1 = math.sqrt(math.pow(carbon1[2]-carbon2[2],2))
    carbon1[2] = carbon1[2]-repeat
    distance2 = math.sqrt(math.pow(carbon1[2]-carbon2[2],2))
    carbon1[2] = carbon1[2] + repeat
    carbon2[2] = carbon2[2] - repeat
    distance3 = math.sqrt(math.pow(carbon1[2]-carbon2[2],2))
    mindistance =  min(distance1,distance2,distance3)
    file.write('\t'+str(mindistance/2)+'\n')
    
    


def setup(folder):
    try:
        os.chdir(folder)
        run()
    except:
        file.write(folder+'\t'+"ERROR\n")
        os.chdir(maindir)

for label in elementsList.keys():
    folder = label + ":"
    setup(folder)
    os.chdir(maindir)
    folder = label + ":" + label
    setup(folder)
    os.chdir(maindir)

file.close
