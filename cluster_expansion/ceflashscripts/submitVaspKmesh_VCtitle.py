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
from kmeshroutines import nstrip, readposcar,create_poscar,readfile,writefile
import dynamicPacking6

def getVCmesh(dir,method,targetNmesh,meshtype):
    lastDir = os.getcwd()   
    meshc = dynamicPacking6.dynamicPack() #instance
    [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = readposcar('POSCAR',dir)
#         create_poscar('POSCAR',descriptor, scale, latticevecs, natoms, postype, positions, path) #just to remove the scale problem
    os.chdir(dir)
    totatoms = sum(natoms)
    atype = 1
    aTypes = []
    for natom in natoms:
        for i in range(natom):
            aTypes.append(atype)
        atype += 1
    aTypes = array(aTypes)
    statusOK,nops = meshc.pack(latticevecs,reciplatt,totatoms,aTypes,postype,transpose(positions),targetNmesh,meshtype,dir,method)
    os.chdir(lastDir)
    return statusOK,nops
    
    
def writejobfile(path,n,type):
    '''read from a template in maindir, and put dir in job name'''
    file1 = open(path +'vaspjob','r')
    jobfile = file1.readlines()
    file1.close
    for i in range(len(jobfile)):
        jobfile[i]=jobfile[i].replace('myjob', dir+'_%i_%s' %(n,type))
    file2 = open(path+'/'+'vaspjob','w')
    file2.writelines(jobfile) 
    file2.close()
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
        os.system ('cp {}/vaspjob {}'.format(vaspinputdir,structDir))  
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
    os.system('cp {}/vaspjob {}'.format(path,newdir))
    writejobfile(newdir,n,type)
    modIncar(newdir)
    return newdir
   

################# script #######################

# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestCUB'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestBCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestFCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestNoMoveFCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestComm/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1DP0.5offset/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/test'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/semicond'
# poscarsDir = '/fslhome/bch/cluster_expansion/vcmesh/semicondPoscars/POSCARS'
# vaspinputdir = '/fslhome/bch/cluster_expansion/vcmesh/semicondPoscars/vaspinput/'
maindir = '/fslhome/bch/cluster_expansion/vcmesh/the99sym'
poscarsDir = '/fslhome/bch/cluster_expansion/vcmesh/poscars99/POSCARS'
vaspinputdir = '/fslhome/bch/cluster_expansion/vcmesh/poscars99/vaspinput'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1059DP/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrBCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrFCC/'
type = 'fcc' 
testfile = 'POSCAR'
method = 0
        #0: exact: use vertices of mesh voronoi cell that are closest/farthest 
        #         from the IBZ center origin to check if the point's volume is cut. 
        #         Cut the VC to determine the volume contribution  
reallatt = zeros((3,3))
# createdirs(poscarsDir,maindir,vaspinputdir)
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
# toRun = []
nopsTot = 0
nstruct = 0
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
            os.chdir(currdir)
#            scale = sum(array(float(poscar[1])))
#            N = rint(Nkppra/sum(array(poscar[5].split(),dtype=int16))).astype(int) # number of kpts desired
#            reallatt[0,:] = array(poscar[2].split())
#            reallatt[1,:] = array(poscar[3].split())
#            reallatt[2,:] = array(poscar[4].split())
#            reallatt = scale*reallatt.astype(float)        
#            reciplatt = 2*pi*transpose(linalg.inv(reallatt))

            os.system('rm slurm*')
            subprocess.call(['rm', 'vasp.out'])
            subprocess.call(['rm', 'OUTCAR'])          
#            subprocess.call(['cp','POSCAR.orig','POSCAR'])
#            subprocess.call(['sbatch', 'vaspjob'])

            for n in range(2,25,1):#23
                print 
                print '==============================================' 
                print 'Base {} in submitVasp (target = n^3)'.format(n)
                print '==============================================' 
                print
                newdir = createRunDir(currdir,n,type) 
                statusOK,nops = getVCmesh(newdir,method,n**3,type)
                nopsTot += nops
                nstruct += 1
                if not statusOK: #no points or too many in IBZ
                    print 'Zero or too many points in IBZ...skip this n'
                    os.system('rm -r {}'.format(newdir))
                else:
                    os.chdir(newdir)
                    subprocess.call(['sbatch', 'vaspjob'])
                    
#                     toRun.append(newdir)
                    
#                 getVCmesh(newdir,method,n,type)
#         sys.exit('stop')
#         newdirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]) 
# for newdir in toRun:
#     os.chdir(newdir)
#     subprocess.call(['sbatch', 'vaspjob'])
#     os.chdir(currdir)
nsymmAvg = nopsTot/float(nstruct)
print 'Average number of symmetry operations:',nsymmAvg, 'for',nstruct,'structures'
os.chdir(maindir)                 
print 'Done'