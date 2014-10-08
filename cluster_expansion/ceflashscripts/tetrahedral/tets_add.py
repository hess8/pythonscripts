
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

class tetrahedrons:
    def __init__(self):
        self.ntet = 0
        self.nbands = 0       
        self.wgt = [] #numpy array
        self.vecs = [] #numpy array
        self.idk = [] ##numpy array; connection table of corners to k's in irreducible BZ. 
        self.vol = [] #volume in units of reciprocal lattice vectors
        
    def remove(self,it):
        self.ntet += -1
        delete(self.wgt,s_[it,:])
        delete(self.vecs,s_[:,:,it])
        delete(self.idk,s_[:,it])
    
    def add(self,tet1):
        print shape(self.idk)
        print shape(tet1.idk)
        tet1.wgt = reshape(tet1.wgt,(1,shape(tet1.wgt)[0]))
        self.wgt = append(self.wgt,tet1.wgt, axis = 0)
        tet1.vecs = reshape(tet1.vecs,(shape(tet1.vecs)[0],shape(tet1.vecs)[1],1))
        append(self.vecs,tet1.vecs, axis = 2)
        append(self.idk,tet1.idk, axis = 1)
    
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
        self.vol = zeros(self.ntet,dtype = float64)
        for it in range(self.ntet):
            for ic in range(4):
                self.vecs[:,ic,it] = lines[5*it+ic+1].split()
            va = self.vecs[:,1,it]-self.vecs[:,0,it]
            vb = self.vecs[:,2,it]-self.vecs[:,0,it]
            vc = self.vecs[:,3,it]-self.vecs[:,0,it]
            self.vol[it] = abs(dot(va,cross(vb,vc))/2.0)
    
    def rd_ibzkpt(self,dir): #only extracts the tets ids
        self.idk = zeros((4,self.ntet),dtype = int)
        lines = nstrip(readfile(dir+'IBZKPT'))
        nk = int(lines[1].split()[0])
        ntet2 = int(lines[4+nk].split()[0])
        if ntet2 != self.ntet: print 'Warning! NTET in IBZKPTS differs from this program''s'
        for it in range(self.ntet):
            self.idk[:,it] = lines[5+nk+it].split()[1:]  #skip the weights           
             
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
    
    def rd_ibzkpt(self,dir): #only extracts the kpoint vectors
        self.vecs = zeros((3,self.nk),dtype = float64)
        lines = nstrip(readfile(dir+'IBZKPT'))
        nk2 = int(lines[1].split()[0])
        if nk2 != self.nk: print 'Warning! Nkpts in IBZKPTS differs from this program''s'
        for ik in range(self.nk):
            self.vecs[:,ik] = lines[ik+3].split()[:3] #don't read weights'
                 
################# script #######################
dir = '/fslhome/bch/cluster_expansion/alal/test/f1_2/'
tet = tetrahedrons()
kp = kpoints()

set_printoptions(precision=15)

os.chdir(dir)
tet.get_tetweights(dir)
print 'tet weights'
print tet.ntet,tet.wgt
kp.get_eigen(dir)
print 'Energies'
print kp.ener
tet.get_tetvecs(dir)
print 'tet vecs'
print tet.vecs
print 'tet volumes'
print tet.vol
tet.rd_ibzkpt(dir)
print 'tet ids to kpoints'
print tet.idk
kp.rd_ibzkpt(dir)
print 'kpoint vecs'
print kp.vecs
print 'nk',kp.nk
print tet.idk[:,1]

#read in connection table

#find partially occupied tets

tet_partial = []
for it in range(tet.ntet):
    for ib in range (1,tet.nbands): #testing weights against filled first band
        if 0.0 < tet.wgt[it,ib] < tet.wgt[it,0]:
#            print it,ib,tet.wgt[it,ib]
            tet_partial.append(it)            
            break #only takes one band to identify it  
#print tet_partial        
#newk = kpoints()
#newk.nk = length(tet_partial)
#newt = tetrahedrons()
#newt.ntet = 4*length(tet_partial)
#newt.vecs = zeros(3,newt.ntet)
newk = kpoints() #just one member, so we can add it
newt = tetrahedrons() #just one member, so we can add it
#newt.wgt = zeros((1,kp.nbands),dtype = float64)
triads = [[0,1,2],[1,2,3],[2,3,0],[3,0,1]] #four faces of old tetrahedron
for it in tet_partial:   
    #add 4 new tets to end, created by adding tets connecting center and other tets
#    print tet.vecs[:,:,it]
    newk.vecs = zeros(3,dtype = float64)
    newt.vecs = zeros((3,4),dtype = float64)
    for ic in range(4): newk.vecs[:] += tet.vecs[:,ic,it]/4.000000000000000  #new center
    for triad in triads:
        newt.vecs[:,0] = newk.vecs[:]
#        print newt.vecs[:,0]

#        print triad
        for ic in range(1,4): 
#            print triad[ic-1]
            newt.vecs[:,ic] = tet.vecs[:,triad[ic-1],it] #three corners for new tet           
#            print newt.vecs[:,ic]
        va = newt.vecs[:,1]-newt.vecs[:,0]
        vb = newt.vecs[:,2]-newt.vecs[:,0]
        vc = newt.vecs[:,3]-newt.vecs[:,0]
#        print 'va, vb, vc',va, vb, vc
        newt.vol = abs(dot(va,cross(vb,vc))/2.0)
#        print 'newt.vol',newt.vol
        newt.wgt = tet.wgt[it,:]*newt.vol/tet.vol[it]
        print 'newt.wgt',newt.wgt
        newt.idk = [kp.nk+1].append([tet.idk[ic] for ic in triad])    
        #still need ids XXXX       
        tet.add(newt)
    #still need ids
#    kp.add(newk)
print 'tet weights'
print tet.ntet,tet.wgt
print 'Energies'
print kp.ener
print 'tet vecs'
print tet.vecs
print 'tet volumes'
print tet.vol
print 'tet ids to kpoints'
print tet.idk
print 'kpoint vecs'
print kp.vecs
    
    
    #remove old tet from structure
#    tet.remove(it-i)  
    #add center to kpoints list
    #update connection table in terms of kpoints list. 
    #      
print 'Done'

