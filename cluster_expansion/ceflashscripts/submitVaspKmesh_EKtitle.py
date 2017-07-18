#!/usr/bin/python

'''Equivalent k point method.  For each POSCAR in poscarsDir in maindir: makes a directory, copies POSCAR, makes a POTCAR from 
the element name in the POSCAR title (double the first atom).  copies vaspinput, creates kpoints file with correct mesh 
(just relative to the recip lattice), reads a jobfile from the vaspinput,
writes the structure tag etc to the job name, and submits a vasp job.  
'''
    
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint
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
    
    
def writekpts_vasp_n(path,n,type):
    '''Write mesh vectors to kpoints file, using integer division for cubic and fcc meshes
    (equivalent k method)'''   
    file1 = open(path +'KPOINTS','w')
    kpointsfile = []
    if type == 'cubic': kpointsfile.append('%i kpoints for cubic integer n=%i\n' %(n**3,n))
    if type == 'fcc': kpointsfile.append('%i kpoints for fcc integer n=%i\n' %(4*n**3,n))
    kpointsfile.append('0 \n')   
    kpointsfile.append('Cart \n')
    print type
    if type == 'cubic':
        b = 1.0/n
        kpointsfile.append('%12.8f 0.0 0.0\n' % b)
        kpointsfile.append('0.0 %12.8f 0.0\n' % b)
        kpointsfile.append('0.0 0.0 %12.8f\n' %  b)
    elif type == 'fcc':
        b = 1.0/n/2.0
        kpointsfile.append('0.0 %12.8f %12.8f\n' % (b,b))
        kpointsfile.append('%12.8f 0.0 %12.8f\n' % (b,b))
        kpointsfile.append('%12.8f %12.8f 0.0\n' % (b,b))
    elif type == 'bcc':
        b = 1.0/n/2.0
        kpointsfile.append('0.0 %12.8f %12.8f\n' % (b,b))
#         kpointsfile.append('%12.8f 0.0 %12.8f\n' % (b,b))
#         kpointsfile.append('%12.8f %12.8f 0.0\n' % (b,b))
    else:
        sys.exit('Stop: type not found in writekpts_vasp_n')
    kpointsfile.append('0.5 0.5 0.5\n' ) #shift
#     kpointsfile.append('0.0 0.0 0.0\n' ) #shift
    file1.writelines(kpointsfile) 
    file1.close()
    return 

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

def createRunDir(path,n,type):
    newdir = path + '%s_%i/' % (type,n)
    if not os.path.isdir(newdir):
        os.system('mkdir %s' % newdir)
    os.system('cp {}/INCAR {}'.format(path,newdir))
    os.system('cp {}/POSCAR {}'.format(path,newdir))
    os.system('cp {}/POTCAR {}'.format(path,newdir))
    os.system('cp {}/vaspjob {}'.format(path,newdir))
    writekpts_vasp_n(newdir,n,type)
    writejobfile(newdir,n,type)  
    return newdir
   

################# script #######################

# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestCUB'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestBCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestFCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestNoMoveFCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestComm/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1DP0.5offset/'
maindir = '/fslhome/bch/cluster_expansion/ekmesh/ekPure'
poscarsDir = '/fslhome/bch/cluster_expansion/ekmesh/pureCubicPoscars/POSCARS'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1059DP/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrBCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrFCC/'
type = 'fcc' 
testfile = 'POSCAR'
vaspinputdir = '/fslhome/bch/cluster_expansion/ekmesh/pureCubicPoscars/vaspinput/'
        #0: exact: use vertices of mesh voronoi cell that are closest/farthest 
        #         from the IBZ center origin to check if the point's volume is cut. 
        #         Cut the VC to determine the volume contribution  
reallatt = zeros((3,3))
createdirs(poscarsDir,maindir,vaspinputdir)
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])

print 'Equivalent k method'
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

            for n in range(2,23,1):#23
                print 
                print '==============================================' 
                print 'Base {} in equiv k submitVasp (target = n^3)'.format(n)
                print '==============================================' 
                print
                newdir = createRunDir(currdir,n,type) 
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
os.chdir(maindir)                 
print 'Done'