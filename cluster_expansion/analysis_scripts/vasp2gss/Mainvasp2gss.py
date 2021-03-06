'''
Starts with completed vasp files, fits and performs ground state search.  No loop
'''

import os, subprocess,sys
from random import seed
from numpy import zeros,sqrt,std,amax,amin
from copy import deepcopy
sys.path.append('/bluehome2/bch/pythonscripts/Erik_scripts/Graphener-d77cc/graphener/') 
import MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo


def contains(struct, alist):
    for i in xrange(len(alist)):
        if str(struct) == str(alist[i]):
            return True
    
    return False

def equals(alist, blist):
    clist = deepcopy(alist)
    dlist = deepcopy(blist)
    while len(clist) > 0:
        if len(dlist) == 0:
            return False
        if not contains(clist[0], dlist):
            return False
        else:
            clist.remove(clist[0])
            dlist.remove(clist[0])
    
    if len(dlist) > 0:
        return False
    else:
        return True
    
def getDiffE(priority, elast, atoms):
    '''Finds the L1 norm energy change between iterations, weighted by priority'''
    ediff = priority[:,:]['energy'] - elast
    ediffL1 = zeros(len(atoms))
    for iatom in range(len(atoms)):
        priorsum = sum(priority[iatom,:]['prior'])
        ediff[iatom,:] = abs(ediff[iatom,:]) * priority[iatom,:]['prior']/priorsum
        ediffL1[iatom] = sqrt(sum(ediff[iatom,:]))
    return ediffL1

def readSettingsFile():
    currDir = os.getcwd()
    infile = open(currDir + '/needed_files/settings.in', 'r')
    inlines = []
    for line in infile:
        firstPart = line.strip().split()[0]
        firstChar = list(firstPart)[0]
        if firstChar != '#':
            inlines.append(line.strip())
    infile.close()
    
    atoms = []
    volRange = []
    clusterNums = []
    trainStructs = 0
    fitStructs = 0
    fitSubsets = 0
    plotTitle = "title"
    xlabel = "xlabel"
    ylabel = "ylabel"
    
    for line in inlines:
        if line.split()[0] == 'ATOMS:':
            i = 1
            adatoms = line.split()[1:]
            for adatom in adatoms:
                if adatom == '#':
                    break
                else:
                    atoms.append(adatom)
            
        elif line.split()[0] == 'VOL_RANGE:':
            low = int(line.split()[1])
            high = int(line.split()[2])
            volRange = [low, high]
            
        elif line.split()[0] == 'CLUSTER_NUMS:':
            parts = line.split()
            for i in xrange(1, 11):
                clusterNums.append(int(parts[i]))
        
        elif line.split()[0] == 'TRAINING_STRUCTS:':
            trainStructs = int(line.split()[1])
          
        elif line.split()[0] == 'FITTING_STRUCTS:':
            fitStructs = int(line.split()[1])
            
        elif line.split()[0] == 'STRUCT_SUBSETS:':
            fitSubsets = int(line.split()[1])
        
        elif line.split()[0] == 'PLOT_TITLE:':
            plotTitle = line.split('\'')[1]
        
        elif line.split()[0] == 'XLAB:':
            xlabel = line.split('\'')[1]
        
        elif line.split()[0] == 'YLAB:':
            ylabel = line.split('\'')[1]
    
    return [atoms, volRange, clusterNums, trainStructs, fitStructs, fitSubsets, plotTitle, xlabel, ylabel]
         
if __name__ == '__main__':
    seed()
    
#    dir = '/fslhome/bch/cluster_expansion/graphene/testLi/'
##    dir = '/fslhome/bch/cluster_expansion/graphene/Li_Al/'
#    os.chdir(dir)

#    maindir = os.getcwd()
    [atoms, volRange, clusterNums, trainingStructs, fitStructs, fitSubsets, plotTitle, xlabel, ylabel] = readSettingsFile()
    uncleOutput = open('uncle_output.txt','w') # All output from UNCLE will be written to this file.    
    ntot = 17010 #replace with enumerator method #bch
    elast = zeros((ntot,len(atoms)),dtype=float) #energies of last iteration, sorted by structure name
    newStructs = []
    gssStructs = []
    lowestStructsFile = open('lowest_vasp.txt','w')
    lowestGssFile = open('lowest_gss.txt','w')
    failedFile = open('failed_vasp.txt','w')
    
    subprocess.call(['echo','\n========================================================'])
    subprocess.call(['echo','Creating structures.in, .holdout, fitting and gss'])
    subprocess.call(['echo','========================================================\n'])
    iteration = 1 #only one 'iteration'


    # Create structures.in and structures.holdout files for each atom.
    uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atoms)
    uncleFileMaker.makeUncleFiles()
    
    # Get all the structs that have been through VASP calculations for each atom. These
    # should be sorted by formation energy during the work done by makeUncleFiles()
    [vStructsFinished, newlyFailed] = uncleFileMaker.getStructureList()
    structuresInLengths = uncleFileMaker.getStructuresInLengths() 
    
    # Perform a fit to the VASP data in structures.in for each atom.
    fitter = Fitter.Fitter(atoms, fitStructs, fitSubsets, structuresInLengths, uncleOutput)
    fitter.makeFitDirectories()
    fitter.fitVASPData(iteration)

    # Perform a ground state search on the fit for each atom.    
    gss = GSS.GSS(atoms, volRange, plotTitle, xlabel, ylabel, uncleOutput)
    gss.makeGSSDirectories()
    gss.performGroundStateSearch(iteration)
    gss.makePlots(iteration)
    
    gss.getAllGSSStructures(iteration, newlyFailed)
    priority = gss.getGssInfo(iteration)
    diffe = getDiffE(priority,elast,atoms)
    print 'diffe',diffe
    elast = priority[:,:]['energy']
    
    uncleOutput.close()
    lowestStructsFile.close()
    lowestGssFile.close()
    failedFile.close()
    print  'Done'
    # Should do some analysis after the loop has finished as well.
    