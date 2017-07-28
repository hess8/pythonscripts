#!/usr/bin/python

'''Dynamic packing method. For each POSCAR in poscarsDir in maindir: makes a directory, copies POSCAR, makes a POTCAR from 
the element name in the POSCAR title (double the first atom).  copies vaspinput, creates kpoints file with correct mesh 
(just relative to the recip lattice), reads a jobfile from the vaspinput,
writes the structure tag etc to the job name, and submits a vasp job.
'''
    
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint
from copy import copy, deepcopy
from numpy.linalg import norm
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/symmetry_k_mesh_search')
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar,create_poscar,readfile,writefile,waitMaxJobs
import dynamicPacking6

def writeJob(path,ntarget,type):
    """ Creates a standard job file for submitting a VASP job on the supercomputer. 
    The job file calls python for mesh definition.
    Writes name and and ntarget."""  
    dir = os.getcwd()
    runFolder = dir.split('/')[-1] 
    jobName = '{}.{}'.format(path[-12:],runFolder)
    jobFile = open('{}/job'.format(path),'w')   
    jobFile.write("#!/bin/bash\n\n")
    jobFile.write('#SBATCH --time=24:30:02\n')
    jobFile.write("#SBATCH --ntasks=4\n")
    jobFile.write("#SBATCH --mem-per-cpu=2G\n")
    jobFile.write("#SBATCH --job-name={}\n".format(jobName)) 
    jobFile.write('module unload mpi\n')
    jobFile.write('module load mpi/openmpi-1.6.5_intel-13.0.1\n')
    jobFile.write('module unload mkl\n')
    jobFile.write('module load mkl/11/2\n')
    jobFile.write('module unload python\n')
    jobFile.write('module load python/2/7\n')
    
    jobFile.write('python ~/pythonscripts/cluster_expansion/ceflashscripts/poscar2meshVC.py {} {} > out\n'.format(ntarget,type))

    jobFile.write('if [ -e "OK" ]\n')
    jobFile.write('then\n')
    jobFile.write('mpiexec ./vasp.x > vasp.out\n')
    jobFile.write('fi\n')
    jobFile.close()
    return

def createdirs(poscarsDir,maindir,vaspinputdir):
    '''makes dir in maindir for each structure in poscarsDir'''
    potcarDir = "/fslhome/bch/vaspfiles/src/potpaw_PBE"
    for file in os.listdir(poscarsDir):
        info = file.split('_')
        struct = '_'.join(info[1:3])
        structDir = '{}/{}'.format(maindir,struct)
        atom = info[1]
        if not os.path.isdir(structDir):
            os.system('mkdir {}'.format(structDir)) #structure is in 
        os.system('cp {}/{} {}/POSCAR'.format(poscarsDir,file,structDir))
        print 'atom',atom
        try:
            potcar = readfile('{}/{}/POTCAR'.format(potcarDir,atom))
        except:
            try:
                potcar = readfile('{}/{}_pv/POTCAR'.format(potcarDir,atom))
            except:
                potcar = readfile('{}/{}_sv/POTCAR'.format(potcarDir,atom))
        potcar += potcar #multiple atoms in POSCAR have only one POTCAR, so they are really all pure cases    
        writefile(potcar,'{}/POTCAR'.format(structDir))
        os.system ('cp {}/INCAR {}'.format(vaspinputdir,structDir))  
        modIncar(structDir)  
        os.system ('cp -P {}/vasp.x {}'.format(vaspinputdir,structDir)) 
    return

def modIncar(newdir):
    lines = readfile('{}/INCAR'.format(newdir))
    plines = readfile('{}/POTCAR'.format(newdir))
    for line in plines:
        if 'ENMAX' in line:
            enmax = 1.4*float(line.split()[2].replace(';',''))
            break
    for i,line in enumerate(deepcopy(lines)):
        if 'ENMAX' in line:
            lines[i] = 'ENMAX = {}\n'.format(enmax)
            break
    writefile(lines,'{}/INCAR'.format(newdir))     
    return

def createRunDir(path,n,type):
    newdir = path + '%s_%i/' % (type,n)
    if not os.path.isdir(newdir):
        os.system('mkdir %s' % newdir)
    os.system('cp {}/INCAR {}'.format(path,newdir))
    os.system('cp {}/POSCAR {}'.format(path,newdir))
    os.system('cp {}/POTCAR {}'.format(path,newdir))
    os.system ('cp -P {}/vasp.x {}'.format(path,newdir))
    writeJob(newdir,n**3,type)
    return newdir
   

################# script #######################

# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestCUB'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestBCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestFCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestNoMoveFCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestComm/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1DP0.5offset/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/test'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/scond_vc'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/scondvr_wc30'
# poscarsDir = '/fslhome/bch/cluster_expansion/vcmesh/semicond/0_info/POSCARS'
# vaspinputdir = '/fslhome/bch/cluster_expansion/vcmesh/semicond/0_info/vaspinput'
maindir = '/fslhome/bch/cluster_expansion/vcmesh/scond_fccOut'
# poscarsDir = '/fslhome/bch/cluster_expansion/vcmesh/the99sym_newMethod/0-info/POSCARS'
# vaspinputdir = '/fslhome/bch/cluster_expansion/vcmesh/the99sym_newMethod/0-info/vaspinput'

# maindir = '/fslhome/bch/cluster_expansion/vcmesh/vr_wc20'
poscarsDir = '{}/0-info/POSCARS'.format(maindir)
vaspinputdir = '{}/0-info/vaspinput'.format(maindir)
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/test2'
# poscarsDir = '/fslhome/bch/cluster_expansion/vcmesh/test2/info/POSCARS'
# vaspinputdir = '/fslhome/bch/cluster_expansion/vcmesh/test2/info/vaspinput'

# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1059DP/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrBCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrFCC/'
type = 'bcc' 
testfile = 'POSCAR' 
reallatt = zeros((3,3))
createdirs(poscarsDir,maindir,vaspinputdir)
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and 'info' not in d])
print 'Dynamic packing method'
for dir in dirs:
    os.chdir(maindir)
    if testfile in os.listdir(dir): 
        print
        currdir = maindir + '/'+ dir+'/'
        print "*********************************************************************************************************"
        print dir + "*********************************************************************************************************"
        print "*********************************************************************************************************"
        file1 = open(currdir+testfile,'r')
        poscar = file1.readlines()
        file1.close()
        if len(poscar) > 0:
#             os.chdir(currdir)
#            scale = sum(array(float(poscar[1])))
#            N = rint(Nkppra/sum(array(poscar[5].split(),dtype=int16))).astype(int) # number of kpts desired
#            reallatt[0,:] = array(poscar[2].split())
#            reallatt[1,:] = array(poscar[3].split())
#            reallatt[2,:] = array(poscar[4].split())
#            reallatt = scale*reallatt.astype(float)        
#            reciplatt = 2*pi*transpose(linalg.inv(reallatt))

#             os.system('rm slurm*')
#             subprocess.call(['rm', 'vasp.out'])
#             subprocess.call(['rm', 'OUTCAR'])          
#            subprocess.call(['cp','POSCAR.orig','POSCAR'])
#            subprocess.call(['sbatch', 'vaspjob'])

            for n in range(2,35,1):#23
                print 
                print '==============================================' 
                print 'Base {} in submitVasp (target = n^3)'.format(n)
                print '==============================================' 
                print
                newdir = createRunDir(currdir,n,type) 
                os.chdir(newdir)
#                 waitMaxJobs()
                subprocess.call(['sbatch', 'job'])
                os.chdir(maindir)
                    
               
print 'Done'