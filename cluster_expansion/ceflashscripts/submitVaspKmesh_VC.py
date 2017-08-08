#!/usr/bin/python

'''For each dir in maindir: copies vaspinput, creates kpoints file with correct mesh 
(just relative to the recip lattice),reads a jobfile from the vaspinput,
writes the structure number etc to the job name, and submits a vasp job.  
Assumes we already have POSCAR from aconvaspPoscar.py
'''
    
import sys,os,subprocess
from pylab import figure,cm,plot,xlabel,ylabel,title,loglog,semilogy,legend,rcParams,style,rc,close
from numpy import zeros,transpose,array,sum,float64,rint
from numpy.linalg import norm
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/symmetry_k_mesh_search')
#import kmeshroutines as km
from kmeshroutines import readposcar,readfile

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
    statusOK = meshc.pack(latticevecs,reciplatt,totatoms,aTypes,postype,transpose(positions),targetNmesh,meshtype,dir,method)
    os.chdir(lastDir)
    return statusOK
    
    
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

def createdir(path,n,type):
    newdir = path + '%s_%i/' % (type,n)
    if not os.path.isdir(newdir):
        os.system('mkdir %s' % newdir)
    os.system ('cp %s* %s' % (vaspinputdir,newdir))
    os.system ('cp %sPOSCAR %s' % (path,newdir))  
    writejobfile(newdir,n,type)  
    return newdir


################# script #######################

# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestCuts'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestBCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cuptTestFCC'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestNoMoveFCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestCuts/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1DP0.5offset/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/the99/'
maindir = '/fslhome/bch/cluster_expansion/vcmesh/test4/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cubicAnalytic/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/12fstrDP_fcc/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/f1059DP/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrBCC/'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/cubicTestRedistrFCC/'
type = 'fcc' 
testfile = 'POSCAR'
vaspinputdir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/vaspinput/'
method = 0.5
        #0: exact: use vertices of mesh voronoi cell that are closest/farthest 
        #         from the IBZ center origin to check if the point's volume is cut. 
        #         Cut the VC to determine the volume contribution  
reallatt = zeros((3,3))
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
# toRun = []
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
            analy=open('analytic','w')
            analy.close()

            for n in range(4,21,1):#23
                print 
                print '==============================================' 
                print 'Target npoints: {}^3'.format(n)
                print '==============================================' 
                print
                newdir = createdir(currdir,n,type) 
                statusOK = getVCmesh(newdir,method,n**3,type)
                if not statusOK: #no points or too many in IBZ
                    print 'Zero or too many points in IBZ...skip this n'
#                     os.system('rm -r {}'.format(newdir))
                else:
                    os.chdir(newdir)
                    subprocess.call(['sbatch', 'vaspjob'])
            os.chdir(currdir)
            lines=readfile('analytic')
            nKs = []
            eners = []
            errs = []
            for line in lines:
                nK = int(line.split()[0])
#                 if nK not in nKs:
                nKs.append(nK)
                eners.append(float(line.split()[1]))
            errs = abs((array(eners) - eners[-1])*1000) + 1e-4
            fig = figure()
            loglog(nKs,errs,marker = 'o',markeredgewidth=0.0,linestyle='None')                             
        
            fig.savefig('analytic_err_vs_n')
            
                
                
os.chdir(maindir)                 
print 'Done'