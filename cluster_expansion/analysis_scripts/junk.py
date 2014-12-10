import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection, symmetryError,\
    latticeType, packingFraction, mink_reduce, lattvec_u,arenormal,\
    unique_anorms, intsymops, create_poscar, searchsphere

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,argmax,nonzero,float64, identity
from numpy import ceil,real,unravel_index, outer, fmod, amin, amax

from scipy.optimize import minimize
from copy import copy,deepcopy
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy.random import randint,random
from itertools import combinations
from pylab import frange
#from ceScriptTools import runningJobs
from pylab import *
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy,vasputil_dosplot


class Fitter:
    """ This class performs the UNCLE fits to the VASP data that has been gathered so far.  It also
        keeps track of the fitting errors, prediction errors, and summaries of cluster expansion 
        terms from iteration to iteration. """

    def __init__(self, atoms, M_fitStructures, N_subsets, vstructsFinished, uncleOutput):
        """ CONSTRUCTOR """  
        self.vstructsFinished = vstructsFinished  
        self.atoms = atoms
        self.M_fitStructures = M_fitStructures
        self.N_subsets = N_subsets
        self.enumFolder = os.getcwd() + '/enum/'
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
        #self.header = "peratom\nnoweights\nposcar\n"
        self.vstructsFinished = vstructsFinished

    def fitVASPData(self, iteration,):
        """ Performs the UNCLE fit to the VASP data. After the fit is done, it adds the iteration
            onto the end of the files we want to keep track of from iteration to iteration. """
        lastDir = os.getcwd()
        
        for iatom, atom in enumerate(self.atoms):

            atomDir = lastDir + '/' + atom
            if os.path.isdir(atomDir):
                subprocess.call(['echo','\nFitting VASP data for ' + atom + '. . .\n'])
                fitsDir = atomDir + '/fits'
                if os.path.isdir(fitsDir):
                    os.chdir(fitsDir)
                    subprocess.call(['uncle', '15']) #not waiting long enough for large cluster numbers

#                        subprocess.call([self.uncleExec, '15'], stdout=self.uncleOut) #not waiting long enough for large cluster numbers
#                        subprocess.call(['mv','fitting_errors.out','fitting_errors_' + str(iteration) + '.out'])
#                        subprocess.call(['mv','prediction_errors.out','prediction_errors_' + str(iteration) + '.out'])
#                        subprocess.call(['mv','J.1.summary.out','J.1.summary_' + str(iteration) + '.out'])
#                        subprocess.call(['cp','structures.in', atomDir]) #so we always have the latest version there, for restarts
#                        subprocess.call(['cp','structures.in', 'structures.in_' + str(iteration)]) #leave the file to be appended to
#                        subprocess.call(['cp','structures.holdout', 'structures.holdout_' + str(iteration)]) #leave the file in case a
#                        os.chdir(lastDir)

os.chdir('/fslhome/bch/cluster_expansion/graphene/tm_row1/')
atoms = ['Ti_sv']
#os.mkdir('junk')
uncleOutput = open('junk/uncle_output.txt','w') # All output from UNCLE will be written to this file. 
fit = Fitter(atoms,0,0,[],uncleOutput)
fit.fitVASPData(0)
for i in range(10):
    print i
uncleOutput.close()