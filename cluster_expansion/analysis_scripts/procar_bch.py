#!/usr/bin/python
'''    
'''
import sys,os,subprocess
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from numpy import zeros,transpose,array,sum,float64,rint,mean,set_printoptions,s_,\
    delete,append,cross,dot
from numpy.linalg import norm
#from numpy.ndarray import flatten
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  getms, writefermi, removezeros

from plotTools import plotxy
from pylab import *
from copy import deepcopy

class procar:
    def __init__(self):
#        self.path = []
        self.nk = 0
        self.kvecs = []
        self.kweights = []
        self.nbands = 0
        self.weights = []
        self.spin = 1 #1: no spin, 2: spin
        self.orbs = []
        self.norbs = 0
        self.nions = 0
        self.ener = []
        self.occ = []
        
    def incar(self,dir):
        print 'reading ISPIN and LORBIT from INCAR'
        lines = nstrip(readfile(dir+'INCAR'))
        for line in lines:
            if ('ISPIN' in line) or ('ispin' in line): self.spin = int(line.split('=')[1])
            if ('LORBIT' in line) or ('lorbit' in line): 
                lorbit = int(line.split('=')[1])
                if lorbit in [0,10]: #no lm decomposition
                    self.orbs = ['s','p','d']
                elif lorbit in [1,2,11,12]:
                    self.orbs = ['s','py','pz','px','dxy','dyz','dz2','dxz','dx2']
                else:
                    sys.exit('Unrecognized LORBIT in INCAR') 
            self.norbs = len(self.orbs)              
        
    def read(self,filepath):
        print 'reading Procar',
        lines = nstrip(readfile(filepath))
        self.nk = int(lines[1].split(':')[1].split()[0])
        self.nbands = int(lines[1].split(':')[2].split()[0])
        self.nions = int(lines[1].split(':')[3])
        #initialize
        self.kvecs = zeros((3,self.nk),dtype = float)
        self.kweights = zeros((self.nk),dtype = float)
        self.ener = zeros((self.nk, self.nbands),dtype = float) 
        self.occ = zeros((self.nk, self.nbands),dtype = float) 
        self.weights = zeros((self.nk,self.nbands,len(self.orbs),self.spin),dtype = float)
        
        #read weights
        ipos = 2
        for ispin in range(self.spin):
            ipos += 1; # print 'ipos for spin %i:  %i' % (ispin,ipos)
            for ik in range(self.nk):
                self.kvecs[:,ik] = lines[ipos].split()[3:6]
                self.kweights[ik] = float(lines[ipos].split()[8])
                ipos += 2
                for ib in range(self.nbands):
                    self.ener[ik,ib] = float(lines[ipos].split()[4])
                    self.occ[ik,ib] = float(lines[ipos].split()[7])
                    ipos += 3
                    for i in range(self.nions):
                        if mod(ipos,1e4)==0: print '-', #progress bar
                        tempw = lines[ipos].split()[1:self.norbs+1]
                        self.weights[ik,ib,:,ispin] = tempw
                        ipos += 1
                    ipos += 2
                ipos += 1
        print       
        
class dos:
    def __init__(self):     
        '''dE is width of gaussian added for each eigenvalue'''
    
    
    def calc(self,pcr):
        emin = amin(pcr.ener.flatten()) 
        emax = amax(pcr.ener.flatten()) 
        print emin, emax
#        dosarr = zeros()                                                                 
################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alal/test/f1_2/'
#maindir = '/fslhome/bch/cluster_expansion/graphene/analysis/1220/DOS/'
maindir = '/fslhome/bch/cluster_expansion/graphene/analysis/1220/DOSlorbit11/'
pcr = procar();dos = dos()
pcr.incar(maindir)
pcr.read(maindir+'PROCAR')

#DOS
dos.calc(pcr)

print 'Done'

