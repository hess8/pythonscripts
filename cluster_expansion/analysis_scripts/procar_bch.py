#!/usr/bin/python
'''    
'''
import sys,os,subprocess
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from numpy import zeros,transpose,array,sum,float64,rint,mean,set_printoptions,s_,\
    delete,append,cross,dot
from numpy.linalg import norm
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
                        print ipos
                        tempw = lines[ipos].split()[1:self.norbs+1]
                        print tempw
                        self.weights[ik,ib,:,ispin] = tempw
                        ipos += 1
                    ipos += 2
                ipos += 1
                                                                                             
################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alal/test/f1_2/'
maindir = '/fslhome/bch/cluster_expansion/graphene/analysis/1220/DOS/'
prcr = procar()
prcr.incar(maindir)
prcr.read(maindir+'PROCAR')

#wrtvasp = write_vasp()
#q = query()
#
#set_printoptions(precision=12)
#os.chdir(maindir)
#dir = '1init/'
#tet.get_tetweights(dir)
#
#kp.get_eigen(dir)
#q.ener(kp)
#tet.get_tetvecs(dir)
##q.tetwgt(tet,kp) #prints
#
##print 'tet8', sum(tet.fwgt[8,:]);tet8=sum(tet.fwgt[8,:])
#print 'tet vecs'
#print tet.vecs
#tet.rd_ibzkpt(dir)
#print 'tet ids to kpoints'
#kp.twgt2kwgt(tet) #kpweights from tetweights
#print tet.idk
#print 'tet volume weights'
#print tet.vwgt
#print 'Sum volume weights', sum(tet.vwgt)
#print 'tet volumes'
#print tet.vol
#kp.rd_ibzkpt(dir)
#print 'kpoint vecs'
#print kp.vecs
#print 'nk',kp.nk
#kp.make_ids(tet)
#print 'kpoint connection to tets'
#print kp.idt
#
#sumfwgt1 = q.sum_fermi(tet,kp)
#add_tets(tet,kp)
#
#sumfwgt2 = q.sum_fermi(tet,kp)
#kp.twgt2kwgt(tet) #kpweights from tetweights
##print 'tet weights'
##print tet.ntet,tet.fwgt
##print 'Energies'
##print kp.ener
##print 'tet vecs'
##print tet.vecs
##print 'tet volume weights'
##print tet.vwgt
##print 'tet volumes'
##print tet.vol
##print 'tet ids to kpoints'
##print tet.idk
##print 'kpoint vecs'
##print kp.vecs
##print 'kpoint to tet ids'
##kp.make_ids(tet)
##print kp.idt
##print 'nk',kp.nk
##print 'ntet',tet.ntet
##print 'Sum volume weights', sum(tet.vwgt)
#
##### VASP run with added k's
#dir = '2addk/'
#wrtvasp.kpfile(tet,kp,dir)  
##run vasp !!
#
#print 'Get new eigenvals, weights, but not vecs'
#tet.get_tetweights(dir)
#
#kp.twgt2kwgt(tet) #kpweights from tetweights
#kp.get_eigen(dir)
#q.ener(kp)
#
#sumfwgt3 = q.sum_fermi(tet,kp)
#
#
##print 'tet 8old', tet8
##print 'tet 26-29' ,sum(tet.fwgt[26:30,:])
##print sumfwgt3
#if abs(sumfwgt3-sumfwgt1)>1e-7: print '  SOMEHOW LOSING FERMI WEIGHT IN VASP!'    
#
##### Add more kpoints
###weights are getting below integers:
##
#add_tets(tet,kp)
#tet.volscale = tet.volscale/4.0
#tet.vwgt = tet.vwgt*4.0
##
###### VASP run with added k's
#dir = '3addk/'
#wrtvasp.kpfile(tet,kp,dir)   
#
#
#
#print 'Done'
