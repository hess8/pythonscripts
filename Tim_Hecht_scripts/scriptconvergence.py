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
file = open('outputconvergencestretch.dat','w')
os.chdir('stretch')

def run():
    
    openfile = open('OSZICAR','r')
    print(folder)
    file.write(folder)
    file.write(openfile.readlines()[-2])

def setup(folder):
    os.chdir(folder)
    run()

for label in elementsList:
    folder = label + ":"
    setup(folder)
    os.chdir("..")
    folder = label + "-" + label + ":"
    setup(folder)
    os.chdir("..")
    folder = label + ":" + label
    setup(folder)
    os.chdir("..")
    folder = label + ":-" + label
    setup(folder)
    os.chdir("..")
    folder = label + "-" + label + ":" + label 
    setup(folder)
    os.chdir("..")
    folder = label + "-" + label + ":" + label + "-" + label
    setup(folder)
    os.chdir("..")

file.close
