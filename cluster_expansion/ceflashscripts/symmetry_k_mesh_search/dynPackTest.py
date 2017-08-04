#!/usr/bin/python
'''Test looping structure for finding the dynamicPacking parameters by fastest 
convergence

All parameters: (Implement search on *'s
        *self.power = 6.0
        *self.wallPower = 6.0
        *self.wallfactor = 1.0  #probably needs to be bigger than interfactor by about the average number of nearest neighbors
        *self.wallClose = 0.1 #0.5 #to allow initial points closer to the wall set to less than 1. 
        *self.wallOffset = 0.5 #back off wall forces and energies by a distance that is a fraction of dw. 
        self.interfactor = 1.0        
        self.initFactor = 1.0
        self.df = 1.00 * self.ravg #inter-point force scale distance
        *self.dw = 0.5 * self.df #wall force scale distance
        eps = self.ravg/300
        *self.searchInitFactor = 0.0
        self.initSrch = 'max'
        self.initSrch = None
        nShift = 5
        nTh = 10
        nPh = 20        
        itermax = 100
        gnormTol = 0.001
        minstep = 0.000001
        step = 1.0 #* minstep
        step /= 2
        neighR = 8.0*self.rpacking'''
    
import sys,os,subprocess
from numpy import (zeros,transpose,array,sum,float64,rint,divide,multiply,argmin,
                   argmax,sign)
from numpy import copy as npcopy
# from copy import copy, deepcopy
from numpy.linalg import norm
from random import uniform,seed
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts')
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/symmetry_k_mesh_search')
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts')
#import kmeshroutines as km
from kmeshroutines import (nstrip,readposcar,create_poscar,readfile,writefile,
                           waitMaxJobs,waitJobs)
import dynamicPacking7, analyzeNks

def setParams(maindir):
    paramLabels = ['power','wallPower','wallfactor','wallClose','wallOffset','dw' ],
    params0 =     [  6.0,     6.0,        1.0,          0.5,         0.5,      0.5 ]
    #to run the gradient in parallel, we need a folder for each parameter.
    for i in range(len(params0)):
        pdir = '{}/p{}'.format(maindir,i)
        if os.path.exists(pdir):
            os.system('rm -r {}'.format(pdir))
        os.mkdir(pdir)
    return params0
    
def Nkcost(params,nlims,maindir,poscarsDir,vaspinputdir):
    createdirs(maindir,poscarsDir,vaspinputdir)
    os.chdir(maindir)
    dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and 'info' not in d])
    jobIDs = []
    ns = []     
    for dir in dirs:
        os.chdir(maindir)
        if testfile in os.listdir(dir): 
            currdir = maindir + '/'+ dir+'/'
            file1 = open(currdir+testfile,'r')
            poscar = file1.readlines()
            file1.close()
            if len(poscar) > 0:
                for n in range(nlims[0],nlims[1],nlims[2]):#23
                    ns.append(n)
                    newdir = createRunDir(currdir,n,type,params) 
                    os.chdir(newdir)
                    waitMaxJobs()
                    proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                    jobid = proc.communicate()[0].split()[3]
                    jobIDs.append(jobid)
                    os.chdir(maindir) 
#     subprocess.call(['echo', '\tSubmitted {} jobs, ID range {} , {} for ns {}'.format(len(jobIDs),jobIDs[0],jobIDs[-1],ns)])
    waitJobs(jobIDs)
    costs = analyzeNks.analyze([maindir])
    return costs[0]

def submit(i,jobIDs,params,nlims,maindir,poscarsDir,vaspinputdir):
    pdir = '{}/p{}'.format(maindir,i)
    createdirs(pdir,poscarsDir,vaspinputdir)
    os.chdir(pdir)
    dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and 'info' not in d])
    ns = []     
    for dir in dirs:
        os.chdir(pdir)
        if testfile in os.listdir(dir): 
            currdir = pdir + '/'+ dir+'/'
            file1 = open(currdir+testfile,'r')
            poscar = file1.readlines()
            file1.close()
            if len(poscar) > 0:
                for n in range(nlims[0],nlims[1],nlims[2]):#23
                    ns.append(n)
                    newdir = createRunDir(currdir,n,type,params) 
                    os.chdir(newdir)
                    waitMaxJobs()
                    proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                    jobid = proc.communicate()[0].split()[3]
                    jobIDs.append(jobid)
                    os.chdir(pdir) 
    return jobIDs

def randSteps(step,x,params0,nlims,maindir,poscarsDir,vaspinputdir):
    '''Advance each component of x separately and return costs and x's. These run in parallel.'''
    jobIDs = []
    xs = []
#     for i in range(len(x)):
#         xtemp = npcopy(x)
#         xtemp[i] *= (1 + sign(uniform(-1,1)) * step )#change only one parameter at a time
#         xs.append(xtemp)
    for i in range(len(x)): #change all parameters with random sign by step
        xtemp = npcopy(x)
        for j in range(len(x)):   
            xtemp[j] *= (1 + sign(uniform(-1,1)) * step )#change only one parameter at a time
        xs.append(xtemp)
        jobIDs = submit(i,jobIDs,multiply(xtemp,params0),nlims,maindir,poscarsDir,vaspinputdir)
    waitJobs(jobIDs)
    rcosts = []
    for i in range(len(x)):
        pdir = '{}/p{}'.format(maindir,i)
        costs = analyzeNks.analyze([pdir])
        rcosts.append(costs[0]) 
    return rcosts,xs

def searchParams(params0,maindir,poscarsDir,vaspinputdir,nlims):
    seed()
    '''The vector x is divide(params,params0).  Deal with x here only.
    Use a random hopping search''' 
    xcurr = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    itermax = 100
    maxSinceMin = 5
    origStep = 0.15
    step = origStep  
    nSinceMin = 1
    nKeep = 10
    best = zeros(nKeep,dtype = [('iIter','int32'),('cost','float'),('x','{}float'.format(len(params0)))])
    best['cost'] += 100
    defCost = Nkcost(multiply(xcurr,params0),nlims,maindir,poscarsDir,vaspinputdir)
    print '\tStarting cost {:6.2f}'.format(defCost)
    best[0]['cost'] = defCost
    best[0]['x'] = xcurr
    iIter = 0
    while iIter < itermax:
        costs,xs = randSteps(step,xcurr,params0,nlims,maindir,poscarsDir,vaspinputdir)
        for i, cost in enumerate(costs):
            xsStr = '['
            for xi in xs[i]:
                xsStr += ' {:6.2f}'.format(xi) 
            print 'cost {}: {:6.2f} {}]'.format(i,cost,xsStr)
            if cost < min(best['cost']):
                nSinceMin = 0 
            if cost < max(best['cost']):
                ibmax = argmax(best['cost'])
                best[ibmax]['cost'] = cost
                best[ibmax]['x'] = xs[i]  
        xsStr = '['
        imin = argmin(costs)
        if nSinceMin == 0:
            step = origStep
            xcurr = xs[imin]
            print 'New x with lower cost {:6.2f}:'.format(costs[imin]),xcurr 
            cStr = ''
            for ic in range(nKeep):
                cStr += ' {:6.2f}'.format(best[ic]['cost']) 
            print '\nLowest costs at iter {}: {}\n'.format(iIter,cStr)
        else:
#             print '\nNo lower cost found at iter {}\n'.format(iIter)
            print 'Iteration {}\n'.format(iIter)
        if nSinceMin > maxSinceMin:
            step *= 2
            print 'Increasing step to', step
            nSinceMin = 1
        iIter += 1 
        nSinceMin += 1                  
#     newParams = currparams
    sort(best,order=['cost']) 
    print 'For {} parameters and {} steps'.format(len(xcurr),iIter)
    print '\tStarting cost {:6.2f}'.format(defCost)
    print '\tLowest cost{:6.2f}'.format(best[0]['cost']),best[0]['x']
    if best[0]['cost'] >= defCost:
        print('Did not find a lower cost')
    if nSinceMin == maxSinceMin:
        sys.exit( '\nStop: no progress during the last {} iterations.'.format(maxSinceMin))
    elif iIter == itermax:
        print '\nExceeded maximum number of iterations ({})'.format(itermax)
    cStr = ''
    for ic in range(nKeep):
        cStr += ' {:6.2f}'.format(best[ic]['cost'])  
    print '\nLowest costs at end{}: {}\n'.format(iIter,cStr)

    return 

def writeJob(path,ntarget,type,params):
    """ Creates a standard job file for submitting a VASP job on the supercomputer. 
    The job file calls python for mesh definition.
    Writes name and and ntarget."""  
    paramStr = ''
    for param in params:
       paramStr += ' {:6.2f}'.format(param)
    dir = os.getcwd()
    runFolder = dir.split('/')[-1] 
    jobName = '{}.{}'.format(path[-12:],runFolder)
    jobFile = open('{}/job'.format(path),'w')   
    jobFile.write("#!/bin/bash\n\n")
    jobFile.write('#SBATCH --time=0:10:02\n')
    jobFile.write("#SBATCH --ntasks=4\n")
    jobFile.write("#SBATCH --mem-per-cpu=2G\n")
    jobFile.write("#SBATCH --job-name={}\n".format(jobName)) 
    jobFile.write('module unload mpi\n')
    jobFile.write('module load mpi/openmpi-1.6.5_intel-13.0.1\n')
    jobFile.write('module unload mkl\n')
    jobFile.write('module load mkl/11/2\n')
    jobFile.write('module unload python\n')
    jobFile.write('module load python/2/7\n')
    
    jobFile.write('python ~/pythonscripts/cluster_expansion/ceflashscripts/poscar2mesh7.py {} {} {} > out\n'.format(ntarget,type,paramStr))

    jobFile.write('if [ -e "OK" ]\n')
    jobFile.write('then\n')
    jobFile.write('mpiexec ./vasp.x > vasp.out\n')
    jobFile.write('fi\n')
    jobFile.close()
    return

def createdirs(maindir,poscarsDir,vaspinputdir):
    '''makes dir in maindir for each structure in poscarsDir'''
    potcarDir = "/fslhome/bch/vaspfiles/src/potpaw_PBE"
    for file in os.listdir(poscarsDir):
        info = file.split('_')
        struct = '_'.join(info[1:3])
        structDir = '{}/{}'.format(maindir,struct)
        atom = info[1]
        if os.path.isdir(structDir):
            os.system('rm -r {}'.format(structDir))
        os.system('mkdir {}'.format(structDir)) #structure is in 
        os.system('cp {}/{} {}/POSCAR'.format(poscarsDir,file,structDir))
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
    for i,line in enumerate(npcopy(lines)):
        if 'ENMAX' in line:
            lines[i] = 'ENMAX = {}\n'.format(enmax)
            break
    writefile(lines,'{}/INCAR'.format(newdir))     
    return

def createRunDir(path,n,type,params):
    newdir = path + '%s_%i/' % (type,n)
    if not os.path.isdir(newdir):
        os.system('mkdir %s' % newdir)
    os.system('cp {}/INCAR {}'.format(path,newdir))
    os.system('cp {}/POSCAR {}'.format(path,newdir))
    os.system('cp {}/POTCAR {}'.format(path,newdir))
    os.system('cp -P {}/vasp.x {}'.format(path,newdir))
    writeJob(newdir,n**3,type,params)
    return newdir

################# script #######################
maindir = '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_SiLP'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_lowPrec'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_lowP2'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/mt_AlLP/'
poscarsDir = '{}/0-info/POSCARS'.format(maindir)
vaspinputdir = '{}/0-info/vaspinput'.format(maindir)

nlims = [2,14,1]
# nlims = [8,12,1]

testfile = 'POSCAR' 
reallatt = zeros((3,3))

print 'Varying parameters in dynamicPacking'
print '\t' + maindir
type = 'bcc'
params0 = setParams(maindir)
xbest = searchParams(params0,maindir,poscarsDir,vaspinputdir,nlims)
print 'Final parameters', multiply(xbest,params0)
               
print 'Done'