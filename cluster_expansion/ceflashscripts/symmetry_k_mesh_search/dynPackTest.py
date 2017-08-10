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
    
import sys,os,subprocess,time
from numpy import (zeros,transpose,array,sum,float64,rint,divide,multiply,argmin,
                   argmax,sign,int32)
from numpy import copy as npcopy
# from copy import copy, npcopy
from numpy.linalg import norm
from random import uniform,seed
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts')
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/symmetry_k_mesh_search')
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts')
#import kmeshroutines as km
from kmeshroutines import (nstrip,readposcar,create_poscar,readfile,writefile,
                           waitMaxJobs,waitJobs)
import dynamicPacking7, analyzeNks

#***************************************
#*************  Settings ***************
maindir = os.getcwd()
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_SiLPnarrdw.5'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_lowPrec'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_lowPrand'
# maindir = '/fslhome/bch/cluster_expansion/vcmesh/mt_LPdw.1/'

#maindir default is os.getcwd()
poscarsDir = '{}/0-info/POSCARS'.format(maindir)
vaspinputdir = '{}/0-info/vaspinput'.format(maindir)
type = 'bcc'
# search = 'grad'
# search = 'rand'
search = 'all'
# nlims = [6,9,1]
nlims = [2,23,1]
#***************************************
#***************************************

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
    jobFile.write('#SBATCH --time=00:40:00\n')
    jobFile.write("#SBATCH --ntasks=8\n")
    jobFile.write("#SBATCH --mem-per-cpu=2G\n")
    jobFile.write("#SBATCH --job-name={}\n".format(jobName)) 
#     jobFile.write("#SBATCH --qos=test\n")
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

def setParams(maindir):
#     paramLabels = ['power','wallPower','wallfactor','wallClose','wallOffset','dw' ]
#     params0 =     [  6.0,     6.0,        1.0,          0.5,         0.5,      0.5 ]
    paramLabels = ['power','wallPower','wallfactor','wallClose','wallOffset' ]
    params0 =     [  6.0,     6.0,        1.0,         0.5,      0.5 ]  
    print 'Initial settings'
    print'\t{}'.format(paramLabels)
    print'\t{}'.format(params0)  
    #to run the gradient in parallel, we need a folder for each parameter.
    for i in range(len(params0)):
        pdir = '{}/p{}'.format(maindir,i)
        if os.path.exists(pdir):
#             os.system('rm -r -f {}'.format(pdir))
            output = subprocess.check_output(['rm','-r','-f',pdir])
        os.mkdir(pdir)
    return params0
    
def Nkcost(params,nlims,dir0,poscarsDir,vaspinputdir):
#     createdirs(dir0,poscarsDir,vaspinputdir)
    os.chdir(dir0)
    dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and 'info' not in d])
    jobIDs = []
    ns = []     
    for dir in dirs:
        os.chdir(dir0)
        if 'POSCAR' in os.listdir(dir): 
            currdir = dir0 + '/'+ dir+'/'
            file1 = open(currdir+'POSCAR','r')
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
                    os.chdir(dir0) 
#     subprocess.call(['echo', '\tSubmitted {} jobs, ID range {} , {} for ns {}'.format(len(jobIDs),jobIDs[0],jobIDs[-1],ns)])
    waitJobs(jobIDs)
    [cost,avgnDone] = analyzeNks.analyze([dir0])
    return cost

def submit(i,jobIDs,params,nlims,maindir,poscarsDir,vaspinputdir):
    pdir = '{}/p{}'.format(maindir,i)
    createdirs(pdir,poscarsDir,vaspinputdir)
    os.chdir(pdir)
    dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and 'info' not in d])
    ns = []     
    for dir in dirs:
        os.chdir(pdir)
        if 'POSCAR' in os.listdir(dir): 
            currdir = pdir + '/'+ dir+'/'
            file1 = open(currdir+'POSCAR','r')
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

def grad(step,x,f,params0,nlims,maindir,poscarsDir,vaspinputdir):
    '''Advance each component of x separately and compute each 
    cost to estimate the gradient. These run in parallel.'''
    jobIDs = []
    f1arr = zeros(len(x),dtype=float)
    gradfactor = (1+step)
    x1 = x * gradfactor
    dx = x1-x
    for i in range(len(x)):
        xtemp = npcopy(x)
        xtemp[i] = x1[i] #change only one parameter at a time
        jobIDs = submit(i,jobIDs,multiply(xtemp,params0),nlims,maindir,poscarsDir,vaspinputdir)
    waitJobs(jobIDs)
    gcosts = []
    print
    for i in range(len(x)):
        pdir = '{}/p{}'.format(maindir,i)
        [cost,avgnDone] = analyzeNks.analyze([pdir])
        f1arr[i] = cost
        gcosts.append(cost) 
        print 'gcost {}: {}'.format(i,f1arr[i])
    grad = divide(f1arr-f,dx)
    print 'grad',grad
    minGradcost = min(gcosts)
    if minGradcost < 0.95*f: #use the corresponding x in the search
        imin = argmin(gcosts)
        xnew = npcopy(x)
        xnew[imin] = x1[imin]
        returnList = [xnew,minGradcost]
    else:
        returnList = None
    return grad,returnList

def randSteps(allRands,step,x,params0,nlims,maindir,poscarsDir,vaspinputdir):
    '''Advance each component of x separately and return costs and x's. These run in parallel.'''
    jobIDs = []
    xs = []
#     for i in range(len(x)):
#         xtemp = npcopy(x)
#         xtemp[i] *= (1 + sign(uniform(-1,1)) * step )#change only one parameter at a time
#         xs.append(xtemp)
    for i in range(len(x)): #change all parameters with random sign by step
        for itry in range(1000):
            xtemp = npcopy(x)
            for j in range(len(x)):   
                xtemp[j] *= (1 + sign(uniform(-1,1)) * step )#change only one parameter at a time
            xsStr = '['
            for xi in xtemp:
                xsStr += ' {:6.2f}'.format(xi)
            if xsStr not in allRands:
                break #this one is unique
        xs.append(xtemp)
        jobIDs = submit(i,jobIDs,multiply(xtemp,params0),nlims,maindir,poscarsDir,vaspinputdir)
    waitJobs(jobIDs)
    rcosts = []
    for i in range(len(x)):
        pdir = '{}/p{}'.format(maindir,i)
        [cost,avgnDone] = analyzeNks.analyze([pdir])
        rcosts.append(cost) 
    return rcosts,xs

def submitSet(ir,params,maindir,poscarsDir,vaspinputdir,nlims):
    '''For use with searchParamsAll.  Submits a parameter set run, and returns to check
    on progress later'''
    jobIDs = []
    rdir = '{}/r{}'.format(maindir,ir)
    createdirs(rdir,poscarsDir,vaspinputdir)
    os.chdir(rdir)
    dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and 'info' not in d])
#     print 'dirs',dirs
    ns = []     
    for dir in dirs:
        os.chdir(dir)
        if 'POSCAR' in os.listdir(os.getcwd()): 
            currdir = rdir + '/'+ dir+'/'
            file1 = open(currdir+'POSCAR','r')
            poscar = file1.readlines()
            file1.close()
            if len(poscar) > 0:
                for n in range(nlims[0],nlims[1],nlims[2]):#23
                    ns.append(n)
                    newdir = createRunDir(currdir,n,type,params) 
                    os.chdir(newdir)
#                     waitMaxJobs()
                    proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                    jobid = proc.communicate()[0].split()[3]
                    jobIDs.append(jobid)
                    os.chdir(rdir) 
    return jobIDs 

def searchParamsAll(maindir,poscarsDir,vaspinputdir,nlims):
    '''Set up a N-dimensional grid of parameters, and run through all of it''' 
    paramLabels = ['power','wallPower','wallfactor','wallClose','wallOffset','dw' ]
    print 'Parameters in method'
    print'\t{}'.format(paramLabels)
    print '\twallPower equals power'
#     print '\tdw held at {}'.format(sys.argv[1])

#     params0 =     [ 2.0, 4.0, 6.0, 8.0 ] 
#     params1 =     'duplicate Power for wallPower' 
#     params2 =     [ 0.1, 0.5, 1.0, 2.0]
#     params3 =     [ 0.1, 0.5, 1.0, 2.0]
#     params4 =     [ 0.0, 0.5, 1.0, 2.0]
 
    params0 =     [6.0] 
    params1 =     'duplicate Power for wallPower' 
    params2 =     [ 0.1, 0.2, 0.3, 0.4]
    params3 =     [ 0.1, 0.2, 0.4, 0.5]
    params4 =     [ 0.0, 0.2, 0.5]
    params5 =     [0.5]
    
#     params0 =     [ 4.0,6.0 ] 
#     params1 =     'duplicate Power for wallPower' 
#     params2 =     [ 0.1]
#     params3 =     [ 0.1]
#     params4 =     [ 0.0]
#     params5 =     [0.5]

#     params5 =  [0.1]
#     params5 =  [float(sys.argv[1])]
    nP =len( paramLabels)
    nPsets = len(params0)*len(params2)*len(params3)*len(params4)*len(params5)
    all = zeros(nPsets,dtype = [('cost','float'),('params','{}float'.format(nP))])
    iset = 0
    nRunSlots = min(100,nPsets)#run slots are directories that run a paramSet.  Should be <= nPsets
    slotsJobIDs = [[]]*nRunSlots
    slotsIsets =  zeros(nRunSlots,dtype = int32)
    for i in range(nRunSlots):
        rdir = '{}/r{}'.format(maindir,i)
        if os.path.exists(rdir):
#             os.system('rm -r -f {}'.format(rdir))
            output = subprocess.check_output(['rm','-r','-f',rdir])
        os.mkdir(rdir)
    isetsToStart = range(nPsets)
    isetsDone = []
    for p0 in params0:
        p1 = p0
        for p2 in params2:
            for p3 in params3:
                for p4 in params4:
                    for p5 in params5:
                        params = [p0,p1,p2,p3,p4,p5]                      
                        all[iset]['params'] = params
                        iset += 1
    summary = open('{}/summaryGrid.csv'.format(maindir),'w')
    summary.write('iset,cost,avgDone,power,wallPower,wallfactor,wallClose,wallOffset,dw\n')
    toAnalyze = []
    iwait = 0
    minCost = 100
    bestParams = []
    while len(isetsDone) < nPsets:
        if len(isetsToStart) > 0:
            icurrSet = isetsToStart[0]
            for ir in range(nRunSlots):
                if len(slotsJobIDs[ir]) == 0: #use this slot for next set
                    iwait = 0; print     
                    #start new set
                    params = all[icurrSet]['params']
                    jobIDs = submitSet(ir,params,maindir,poscarsDir,vaspinputdir,nlims)
                    subprocess.call(['echo', '\tFor set {} in slot {}, submitted {} jobs, ID range {} , {}'.format(icurrSet+1,ir,len(jobIDs),jobIDs[0],jobIDs[-1],icurrSet)])
                    slotsJobIDs[ir] = jobIDs
                    isetsToStart.pop(0)
                    toAnalyze.append(ir)
                    slotsIsets[ir] = icurrSet
                    if len(slotsJobIDs[-1]) > 0: #slots have all started work
                        print '\twait', 
                    break #submit one set at a time
        if len(toAnalyze) > 0:
            for ir in range(nRunSlots):
                if len(slotsJobIDs[ir]) == 0: 
                    if ir in toAnalyze: #the slot's previous calc has not been analyzed
                        print
                        setDir = '{}/r{}'.format(maindir,ir)
                        [cost,avgnDone] = analyzeNks.analyze([setDir])
                        ioldSet = slotsIsets[ir]
                        all[ioldSet]['cost'] = cost
                        if cost < minCost: 
                            minCost = cost
                            bestAvgNdone = avgnDone
                            bestParams = all[ioldSet]['params']
                            iminCost =  ioldSet
                            rdir = '{}/r{}'.format(maindir,ir)
                            os.system('mv {}/loglog_e_vs_n.png {}/best_loglog_e_vs_n.png'.format(rdir,maindir))
                            os.system('mv {}/methodErrs.png {}/best_methodErrs.png'.format(rdir,maindir)) 
                            os.system('mv {}/summary.csv {}/best_summary.csv'.format(rdir,maindir))                    
                        print 'cost for set {}: {:6.2f} {}] avg nDone {}'.format(ioldSet,cost,all[ioldSet]['params'],avgnDone)
                        print 'vs. min cost {}: {:6.2f} {}] avg nDone {}'.format(iminCost,minCost,bestParams,bestAvgNdone)
                        ps = all[ioldSet]['params']
                        summary.write('{},{:6.3f},{:6.2f},{},{},{},{},{},{}\n'.format(ioldSet,cost,avgnDone,
                                                    ps[0],ps[1],ps[2],ps[3],ps[4],ps[5]))
                        summary.flush()
                        toAnalyze.remove(ir)
                        isetsDone.append(ioldSet) 
                        break #analyze one set at a tome
    
        #update slotsJobIDs
        output = []
        while len(output) in [0]: #[0,8]:
            devnull = open(os.devnull, 'w')
            proc = subprocess.Popen(['squeue', '-u', 'bch'], stdout=subprocess.PIPE, stderr=devnull)
            output = proc.communicate()[0].split()
        if len(output)>8:
            runningIDs = []
            for item in output:
                if item.isdigit() and len(item) == 8:
                    runningIDs.append(item)
            for ir in range(nRunSlots):
                for id in npcopy(slotsJobIDs[ir]):
                    if id in runningIDs:
                        continue
                    else:
                        slotsJobIDs[ir].remove(id) 
        elif len(output) == 8:
            slotsJobIDs = [[]]*nRunSlots
        if len(slotsJobIDs[-1]) > 0:
               iwait += 1   
               time.sleep(10)
               print ' {}'.format(iwait),
    summary.close()  
    print 'Finished {} sets of parameters'.format(nPsets)
    
def searchParamsRand(params0,maindir,poscarsDir,vaspinputdir,nlims):
    '''The vector x is divide(params,params0).  Deal with x here only.
    Use a random hopping search''' 
    seed()
    xcurr = array([1.0, 1.0, 1.0, 1.0, 1.0])
    itermax = 100
    maxSinceMin = 15
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
    allRands = [] 
    while iIter < itermax:
        costs,xs = randSteps(allRands,step,xcurr,params0,nlims,maindir,poscarsDir,vaspinputdir)
        for i, cost in enumerate(costs):
            xsStr = '['
            for xi in xs[i]:
                xsStr += ' {:6.2f}'.format(xi) 
            allRands.append(xsStr)
            print 'cost {}: {:6.2f} {}]'.format(i,cost,xsStr)
            if cost < min(best['cost']):
                nSinceMin = 0 
            if cost < max(best['cost']):
                ibmax = argmax(best['cost'])
                best[ibmax]['cost'] = cost
                best[ibmax]['x'] = xs[i]  
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

def searchParams(params0,maindir,poscarsDir,vaspinputdir,nlims):
    '''The vector x is divide(params,params0).  Deal with x here only.
    Uses gradient search''' 
#     xbest = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    xbest = array([1.0, 1.0, 1.0, 1.0, 1.0])
    itermax = 100
    gnormTol = 0.001
    minstep = 0.0001
    iIter = 0
    step = 0.1  
    fbest = Nkcost(multiply(xbest,params0),nlims,maindir,poscarsDir,vaspinputdir)
    print '\tDefault parameters cost {:6.2f}'.format(fbest)
    returnList = []
    while not returnList is None:
        if len(returnList)>0:
            print 'Moving to low cost point found in grad routine'
            xbest = returnList[0]
            fbest = returnList[1]
        gr,returnList = grad(max([step,0.01]),xbest,fbest,params0,nlims,maindir,poscarsDir,vaspinputdir)
    gnorm  = norm(gr)
    gnormstart = gnorm
    fstart = fbest
    method = 'steepest'
    print 'Initial cost',fbest,'gnorm',gnorm,xbest 
    atMinStep = False
    while iIter < itermax and gnorm > gnormTol and not atMinStep:
        print iIter, #progress bar   
        lower = False
        while not lower:
            if method == 'steepest':
                xnew = xbest - step*gr
            fnew =  Nkcost(multiply(xnew,params0),nlims,maindir,poscarsDir,vaspinputdir)
            print '\tCost {:6.2f}'.format(fnew),xnew,step
            if fnew < fbest:
                lower = True
                fbest = fnew
                xbest = npcopy(xnew)
                returnList = []
                while not returnList is None:
                    if len(returnList)>0:
                        print 'Moving to low cost point found in grad routine'
                        xbest = returnList[0]
                        fbest = returnList[1] 
                    gr,returnList = grad(max([step,0.01]),npcopy(xbest),fbest,params0,nlims,maindir,poscarsDir,vaspinputdir)
                gnorm  = norm(gr)
                step *= 2
            else:
                step /= 2
                if step < minstep:
                    print 'minimum step reached: {}'.format(step) 
                    atMinStep = True
                    break
        iIter += 1                   
#     newParams = currparams
    print 'For {} parameters and {} steps'.format(len(xbest),iIter)
    print '\tStarting cost',fstart, 'gnorm',gnormstart
    print '\tEnding cost',fnew,'gnorm',gnorm,'step',step#, 'grad', gnew
    if gnorm <= gnormTol:
        print '\nSuccess after {} iterations'.format(iIter)
    elif iIter == itermax:
        print '\nExceeded maximum number of iterations ({}), while gnorm {} is greater than the tolerance {}'.format(itermax,gnorm,gnormTol)
    if fnew >= fstart:
        print('Did not find a lower cost!')


def createdirs(maindir,poscarsDir,vaspinputdir):
    '''makes dir in maindir for each structure in poscarsDir'''
    potcarDir = "/fslhome/bch/vaspfiles/src/potpaw_PBE"
    for file in os.listdir(poscarsDir):
        info = file.split('_')
        struct = '_'.join(info[1:3])
        structDir = '{}/{}'.format(maindir,struct)
        atom = info[1]
        if os.path.isdir(structDir):
#             os.system('rm -r -f {}'.format(structDir))
            output = subprocess.check_output(['rm','-r','-f',structDir])
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

#search = 'rand'
print 'Varying parameters in dynamicPacking'
print '\t' + maindir


if search == 'grad':
    params0 = setParams(maindir)
    searchParams(params0,maindir,poscarsDir,vaspinputdir,nlims)
elif search == 'rand':
    params0 = setParams(maindir)
    searchParamsRand(params0,maindir,poscarsDir,vaspinputdir,nlims)
elif search == 'all':
    searchParamsAll(maindir,poscarsDir,vaspinputdir,nlims)
               
print 'Done'