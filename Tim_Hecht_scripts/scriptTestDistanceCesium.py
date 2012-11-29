#script.py
import os,subprocess
print("Initializing...\n")

#outputFile = open("distanceTest.dat",'w')
potcarDir = "/bluehome/thecht/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/"
mainDir = "/bluehome/jshum33/VASP/MD/Gradualstretching"


ITERATIONS = 10
PROCESSORS = 32

def run():

    proc = subprocess.Popen(['mpirun','-np',str(PROCESSORS),'vasp.mpi'],stdout=subprocess.PIPE)
    while proc.poll() is None:
        output = proc.stdout.readline()
        file.write(output)
        print output,


def setup(distance):
    if checkFolder()==True:
        print("\nPreparing Files...")
        Poscar(distance)
        Incar()
        Potcar()
        Kpoints()
        SymbolicLink()
        run()

i=1

while i < ITERATIONS+1:
    print "Running iteration " + str(i)
    p = subprocess.check_call(['/fslapps/matlab/matlab_7.8/bin/matlab', '-nodisplay' ,  '-nojvm' , '-nosplash' , '-r', 'StretchingProgramDirect75'])
    run()
    p = subprocess.check_call(['mv','POSCAR','POSCAR'+str(i)])
    p = subprocess.check_call(['mv','OUTCAR','OUTCAR'+str(i)])
    p = subprocess.check_call(['mv','OSZICAR','OSZICAR'+str(i)])
    
