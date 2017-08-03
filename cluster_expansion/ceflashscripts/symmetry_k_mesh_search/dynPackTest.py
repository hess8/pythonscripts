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
from random import uniform
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
   
# def cost(x,params0,nlims,maindir,poscarsDir,vaspinputdir):
#     '''Compute current cost'''
#     f = Nkcost(multiply(x,params0),nlims,maindir,poscarsDir,vaspinputdir)
#     return f

# def grad(step,x,f,params0,nlims,maindir,poscarsDir,vaspinputdir):
#     '''Advance each component of x separately and compute each 
#     cost to estimate the gradient. These run in parallel.'''
#     jobIDs = []
#     f1arr = zeros(len(x),dtype=float)
#     gradfactor = (1+step)
#     x1 = x * gradfactor
#     dx = x1-x
#     for i in range(len(x)):
#         xtemp = x
#         xtemp[i] = x1[i] #change only one parameter at a time
#         jobIDs = gsubmit(i,jobIDs,multiply(xtemp,params0),nlims,maindir,poscarsDir,vaspinputdir)
#     waitJobs(jobIDs)
#     gcosts = []
#     for i in range(len(x)):
#         pdir = '{}/p{}'.format(maindir,i)
#         costs = analyzeNks.analyze([pdir])
#         f1arr[i] = costs[0]
#         gcosts.append(costs[0]) 
#         print 'gcost {}: {}'.format(i,f1arr[i])
#     grad = divide(f1arr-f,dx)
#     print '\ngrad',grad
#     minGradcost = min(gcosts)
#     if minGradcost < 0.85*f: #use the corresponding x in the search
#         imin = argmin(gcosts)
#         xnew = deepcopy(x)
#         xnew[i] = x1[i]
#         returnList = [xnew,minGradcost]
#     else:
#         returnList = None
#     return grad,returnList

def randSteps(step,x,params0,nlims,maindir,poscarsDir,vaspinputdir):
    '''Advance each component of x separately and return costs and x's. These run in parallel.'''
    jobIDs = []
    xs = []
    for i in range(len(x)):
        xtemp = npcopy(x)
        xtemp[i] *= (1 + sign(uniform(-1,1)) * step )#change only one parameter at a time
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
    '''The vector x is divide(params,params0).  Deal with x here only.
    Use a random hopping search''' 
    xcurr = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    itermax = 100
    notLowerMax = 20
    step = 0.1  
    nNotLower = 1
    nKeep = 10
    best = zeros(nKeep,dtype = [('iIter','int32'),('cost','float'),('x','{}float'.format(len(params0)))])
    best['cost'] += 100
    defCost = Nkcost(multiply(xcurr,params0),nlims,maindir,poscarsDir,vaspinputdir)
    print '\tStarting cost',defCost
    best[0]['cost'] = defCost
    best[0]['x'] = xcurr
    iIter = 0
    while iIter < itermax and nNotLower < notLowerMax:
        costs,xs = randSteps(step,xcurr,params0,nlims,maindir,poscarsDir,vaspinputdir)
        for i, cost in enumerate(costs):
            xsStr = '['
            for ix in xs[i]:
                xsStr += ' {:6.3f}'.format(ix) 
            print 'cost {}: {:6.3f} xi: {}]'.format(i,cost,xsStr)
            if cost < max(best['cost']):
                nNotLower = 0
                ibmax = argmax(best['cost'])
                best[ibmax]['cost'] = cost
                best[ibmax]['x'] = xs[i]      
        print '\nLowest costs at iter {}:'.format(iIter),best['cost']
        imin = argmin(costs)
        xcurr = xs[imin]
        if nNotLower == 0:
            print 'New x with cost {}:'.format(costs[imin]),xcurr 
        iIter += 1 
        nNotLower += 1                  
#     newParams = currparams
    best.sort('cost')
    print 'For {} parameters and {} steps'.format(len(xcurr),iIter)
    print '\tStarting cost',defCost
    print '\tLowest cost',best[0]['cost'],best[0]['x']
    if best[0]['cost'] >= defCost:
        print('Did not find a lower cost')
    if nNotLower == notLowerMax:
        print '\nStop: no progress during the latest {} iterations.'.format(notLowerMax)
    elif iIter == itermax:
        print '\nExceeded maximum number of iterations ({})'.format(itermax)
    return 

                           
# def searchParams(params0,maindir,poscarsDir,vaspinputdir,nlims):
#     '''The vector x is divide(params,params0).  Deal with x here only.''' 
#     xbest = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
#     itermax = 100
#     gnormTol = 0.001
#     minstep = 0.0001
#     iIter = 0
#     step = 0.1  
#     fbest = cost(deepcopy(xbest),params0,nlims,maindir,poscarsDir,vaspinputdir)
#     print 'Cost at defaults',fbest
#     checkLow = []
#     while not checkLow is None:
#         if len(checkLow)>0:
#             print 'Moving to low cost point found in grad routine'
#             xbest = deepcopy(checkLow[0])
#             fbest = checkLow[1]
#         gr,checkLow = grad(max([step,0.01]),deepcopy(xbest),fbest,params0,nlims,maindir,poscarsDir,vaspinputdir)
#     gnorm  = norm(gr)
#     gnormstart = gnorm
#     fstart = fbest
#     method = 'steepest'
#     print 'Cost before iterations',fbest,'gnorm',gnorm,xbest 
#     atMinStep = False
#     while iIter < itermax and gnorm > gnormTol and not atMinStep:
#         print iIter, #progress bar   
#         lower = False
#         while not lower:
#             if method == 'steepest':
#                 xnew = xbest - step*gr
#             fnew = cost(deepcopy(xnew),params0,nlims,maindir,poscarsDir,vaspinputdir)
#             print 'Cost',fnew,xnew,step
#             if fnew < fbest:
#                 lower = True
#                 fbest = copy(fnew)
#                 xbest = deepcopy(xnew)
#                 checkLow = []
#                 while not checkLow is None:
#                     if len(checkLow)>0:
#                         print 'Moving to low cost point found in grad routine'
#                         xbest = checkLow[0]
#                         fbest = checkLow[1] 
#                     gr,checkLow = grad(max([step,0.01]),deepcopy(xbest),fbest,params0,nlims,maindir,poscarsDir,vaspinputdir)
#                 gnorm  = norm(gr)
#                 step *= 2
#             else:
#                 step /= 2
#                 if step < minstep:
#                     print 'minimum step reached: {}'.format(step) 
#                     atMinStep = True
#                     break
#         iIter += 1                   
# #     newParams = currparams
#     print 'For {} parameters and {} steps'.format(len(xnew),iIter)
#     print '\tStarting cost',fstart, 'gnorm',gnormstart
#     print '\tEnding cost',fnew,'gnorm',gnorm,'step',step#, 'grad', gnew
#     if gnorm <= gnormTol:
#         print '\nSuccess after {} iterations'.format(iIter)
#     elif iIter == itermax:
#         print '\nExceeded maximum number of iterations ({}), while gnorm {} is greater than the tolerance {}'.format(itermax,gnorm,gnormTol)
#     if fnew >= fstart:
#         print('Did not find a lower cost!')
#     return xbest

def writeJob(path,ntarget,type,params):
    """ Creates a standard job file for submitting a VASP job on the supercomputer. 
    The job file calls python for mesh definition.
    Writes name and and ntarget."""  
    paramStr = ''
    for param in params:
       paramStr += ' {:6.3f}'.format(param)
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
maindir = '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_lowPrec'
poscarsDir = '{}/0-info/POSCARS'.format(maindir)
vaspinputdir = '{}/0-info/vaspinput'.format(maindir)

# nlims = [2,14,1]
nlims = [5,21,1]

testfile = 'POSCAR' 
reallatt = zeros((3,3))

print 'Varying parameters in dynamicPacking'
print '\t' + maindir
type = 'bcc'
params0 = setParams(maindir)
xbest = searchParams(params0,maindir,poscarsDir,vaspinputdir,nlims)
print 'Final parameters', multiply(xbest,params0)
               
print 'Done'