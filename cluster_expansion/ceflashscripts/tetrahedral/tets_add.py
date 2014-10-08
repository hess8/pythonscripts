
#!/usr/bin/python
'''    
'''

import sys,os,subprocess
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from numpy import zeros,transpose,array,sum,float64,rint,mean,set_printoptions,s_,\
    delete,append
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  getms, writefermi, removezeros

from plotTools import plotxy
from pylab import *
from copy import deepcopy

class tetrahedrons:
    def __init__(self):
        self.ntet = 0
        self.nbands = 0       
        self.wgt = [] #numpy array
        self.vecs = [] #numpy array
        self.idk = [] ##numpy array; connection table of corners to k's in irreducible BZ. 
    
    def remove(self,it):
        self.ntet += -1
        delete(self.wgt,s_[it,:])
        delete(self.vecs,s_[:,:,it])
        delete(self.idk,s_[:,it])
    
    def add(self,tet1):
        append(self.wgt,tet1.wgt, axis = 0)
        append(self.vecs,tet1.vecs, axis = 2)
        append(self.vecs,tet1.idk, axis = 1)
    
    def get_tetweights(self,dir):
        lines = nstrip(readfile(dir+'tetwgt'))
        self.ntet = int(lines[0].split()[0])
        self.nbands = int(lines[0].split()[1])
        #    nspin = int(lines[0][2]) #haven't added spin yet
        self.wgt = zeros((self.ntet,self.nbands),dtype = float64)
        for it in range(self.ntet):
            for ib in range(self.nbands):
                self.wgt[it,ib] = lines[it+1].split()[ib+1] #First line and first column are not weights                   

    def get_tetvecs(self,dir):
        lines = nstrip(readfile(dir+'tetvecs'))
        self.vecs = zeros((3,4,self.ntet),dtype = float64)
        for it in range(self.ntet):
            for ic in range(4):
                self.vecs[:,ic,it] = lines[5*it+ic+1].split()              
        
class kpoints:
    def __init__(self):
        self.nk = []
        self.nbands = 0   
        self.vecs = []
        self.idtet = []#connection table to tetrahedra it is a corner to.  Each k is part of 4 tets, 
        self.ener = [] 
        
    def add(self,kpt1):
        append(self.vecs,kpt1.vecs, axis = 1)
        append(self.vecs,kpt1.idtet, axis = 1)

    def get_eigen(self,dir):
        lines = nstrip(readfile(dir+'eigenv'))
        self.nk = int(lines[0].split()[0])
        self.nbands = int(lines[0].split()[1])
        self.ener = zeros((self.nk,self.nbands),dtype = float64)
        for ik in range(self.nk):
            for ib in range(self.nbands):
                self.ener[ik,ib] = lines[ik+1].split()[ib+1]             


################# script #######################
dir = '/fslhome/bch/cluster_expansion/alal/test/f1_2/'
tet = tetrahedrons()
kp = kpoints()

set_printoptions(precision=15)

os.chdir(dir)
tet.get_tetweights(dir)
print tet.ntet,tet.wgt
kp.get_eigen(dir)
print kp.ener
tet.get_tetvecs(dir)
print tet.vecs

#read in connection table

#find partially occupied tets

tet_partial = []
for it in range(tet.ntet):
    for ib in range (1,tet.nbands): #testing weights against filled first band
        if 0.0 < tet.wgt[it,ib] < tet.wgt[it,0]:
            print it,ib,tet.wgt[it,ib]
            tet_partial.append(it)            
            break #only takes one band to identify it  
        
newk = kpoints()
newk.nk = length(tet_partial)
newt = tetrahedrons()
newt.ntet = 4*length(tet_partial)
newt.vecs = zeros(3,newt.ntet)
for it_old, i in enumerate(tet_partial):   
    #add 4 new tets to end, created by adding tets connecting center and other tets
     #will just have one member
    newtet.vecs[ = zeros(3)
    for ic in range(4): center += tet.vecs[:,it,ic]/4.0
    
    #remove old tet from structure
    tet.remove(it_old-i)  
    #add center to kpoints list
    #update connection table in terms of kpoints list. 
    #      
print 'Done'

