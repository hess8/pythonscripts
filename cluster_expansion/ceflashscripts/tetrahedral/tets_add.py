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
        self.fwgt = [] #numpy array Fermi weight of each tetrahedron
        self.vwgt = [] #numpy array Volume weight of each tetrahedron        
        self.vecs = [] #numpy array
        self.idk = [] ##numpy array; connection table of corners to k's in irreducible BZ. 
        self.vol = [] #volume in units of reciprocal lattice vectors
        self.volscale = [] #vasp's volume scale
        
    def remove(self,it):
        self.ntet += -1
        self.fwgt = delete(self.fwgt,it,0)
        self.vecs = delete(self.vecs,it,2)
        self.idk = delete(self.idk,it,1)
        self.vol = delete(self.vol,it)
        self.vwgt = delete(self.vwgt,it)
    
    def add(self,tet1):
        arr = reshape(tet1.fwgt,(1,shape(tet1.fwgt)[0]))
        self.fwgt = self.fwgt = append(self.fwgt,arr, axis = 0)
        arr = reshape(tet1.vecs,(shape(tet1.vecs)[0],shape(tet1.vecs)[1],1))
        self.vecs = append(self.vecs,arr,axis = 2)
        arr = reshape(tet1.idk,(shape(tet1.idk)[0],1))
        self.idk = append(self.idk,arr,axis = 1)
        arr = array(tet1.vol)
        self.vol = append(self.vol,arr)
        self.vwgt = append(self.vwgt,tet1.vwgt)        
        self.ntet += 1
    
    def get_tetweights(self,dir):
        lines = nstrip(readfile(dir+'tetwgt'))
        self.ntet = int(lines[0].split()[0])
        self.nbands = int(lines[0].split()[1])
        #    nspin = int(lines[0][2]) #haven't added spin yet
        self.fwgt = zeros((self.ntet,self.nbands),dtype = float64)
        for it in range(self.ntet):
            for ib in range(self.nbands):
                self.fwgt[it,ib] = lines[it+1].split()[ib+1] #First line and first column are not weights                   

    def get_tetvecs(self,dir):
        lines = nstrip(readfile(dir+'tetvecs'))
        self.vecs = zeros((3,4,self.ntet),dtype = float64)
        self.vol = zeros(self.ntet,dtype = float64)
        print shape(self.vecs)
        for it in range(self.ntet):
            for ic in range(4):
                print it, ic,lines[5*it+ic+1]
                self.vecs[:,ic,it] = lines[5*it+ic+1].split()
            va = self.vecs[:,1,it]-self.vecs[:,0,it]
            vb = self.vecs[:,2,it]-self.vecs[:,0,it]
            vc = self.vecs[:,3,it]-self.vecs[:,0,it]
            self.vol[it] = abs(dot(va,cross(vb,vc))/2.0)
    
    def rd_ibzkpt(self,dir): #only extracts the tets ids
        self.idk = zeros((4,self.ntet),dtype = int)
        self.vwgt = zeros(self.ntet,dtype = float64)
        lines = nstrip(readfile(dir+'IBZKPT'))
        nk = int(lines[1].split()[0])
        ntet2 = int(lines[4+nk].split()[0])
        self.volscale = float(lines[4+nk].split()[1])
        print 'VASP volume scale', self.volscale
        if ntet2 != self.ntet: print 'Warning! NTET in IBZKPTS differs from this program''s'
        for it in range(self.ntet):
            self.vwgt[it] = lines[5+nk+it].split()[0]
            self.idk[:,it] = lines[5+nk+it].split()[1:]         
        self.idk = self.idk -1 #(converting from Fortran counting to python)

class kpoints:
    def __init__(self):
        self.nk = []
        self.nbands = 0   
        self.vecs = []
        self.idt = []#connection table to tetrahedra it is a corner to.  Each k is part of 4 tets, 
        self.ener = [] 
        self.degen = []
        self.twgt = [] #weight derived from tetrahedra volume weight, not from symmetry ops on kpoints
        
    def add(self,kpt1):
        self.nk += 1
        arr = reshape(kpt1.vecs,(shape(kpt1.vecs)[0],1))              
        self.vecs = append(self.vecs,arr, axis = 1)
        
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
        if nk2 != self.nk: print nk2,self.nk, 'Warning! Nkpts in IBZKPTS differs from this program''s'
        for ik in range(self.nk):
            self.vecs[:,ik] = lines[ik+3].split()[:3] #don't read weights'

    def make_ids(self,tet):  #creates the connection ids table between this kpoint and the tetrahedrons it is part of
        self.idt = -1 + zeros((tet.ntet,self.nk),dtype=int) #null is designated -1, becuase we use tet 0 as the first one.
        idcount = zeros(self.nk) #each may contribute to 4 tetrahedra. In full BZ is it exactly 4. Here it may be less because of symmetry
        for ik in range(self.nk):
            for it in range(tet.ntet):
                if ik in tet.idk[:,it]: 
                    self.idt[idcount[ik],ik] = it
                    idcount[ik] += 1   
    def twgt2kwgt(self,tet):
        self.twgt = zeros(self.nk,dtype=float)
        for ik in range(self.nk):
            for it in range(tet.ntet):
                for ic in range(4):
                    if ik == tet.idk[ic,it]: 
                        self.twgt[ik] += tet.vwgt[it]/4.0
        for ik in range(self.nk):
            if mod(self.twgt[ik],1) == 0.5:
                self.twgt = self.twgt * 2.0 #so we have integers
                break
        print 'Unscaled k point weights from tets'
        for ik in range(self.nk): print ik, self.twgt[ik]
        print 'Total weight',sum(self.twgt)
                              
class write_vasp:
    def __init__(self):
        ''''''
    def kpfile(self,tet,kp,dir):
        print 'Writing KPOINTS to dir %s' % dir
        f = open(dir+'KPOINTS','w')
        f.write('BCH generated from tets_add.py\n')        
        f.write(str(kp.nk)+'\n')
        f.write('Reciprocal\n')
        for ik in range(kp.nk):
            f.write('%16.12f %16.12f %16.12f %16.12f\n' % (kp.vecs[0,ik],kp.vecs[1,ik],kp.vecs[2,ik],1.0))  #kp.twgt[ik]) )      
        f.write('Tetrahedra\n')
        f.write('%i       %16.12f\n' % (tet.ntet,tet.volscale))
        for it in range(tet.ntet):
            f.write('%f   %i    %i   %i   %i\n'  % (tet.vwgt[it],tet.idk[0,it]+1,tet.idk[1,it]+1,tet.idk[2,it]+1,tet.idk[3,it]+1))
        f.close()

def add_tets(tet,kp):
    #find partially occupied tets
    tet_partial = []
    for it in range(tet.ntet):
        for ib in range (1,tet.nbands): #testing weights against filled first band
            if 0.0 < tet.fwgt[it,ib] < tet.fwgt[it,0]:
    #            print it,ib,tet.fwgt[it,ib]
                tet_partial.append(it)            
                break #only takes one band to identify it  
    newk = kpoints() #just one member, so we can add it
    newt = tetrahedrons() #just one member, so we can add it
    triads = [[0,1,2],[1,2,3],[2,3,0],[3,0,1]] #four faces of old tet with new center become four new tets
    for i,it in enumerate(tet_partial):   
        #add 4 new tets to end, created by adding tets connecting center and other tets
        newt.vecs = zeros((3,4),dtype = float64)
        newk.vecs = zeros(3,dtype = float64)
        newk.idt = zeros(tet.ntet+4,dtype = int) #add 4 tets
        for ic in range(4): newk.vecs[:] += tet.vecs[:,ic,it-i]/4.0 #center is new kpoint
        for triad in triads: #make a new tetrahedron
            newt.vecs[:,0] = newk.vecs[:]
            for ic in range(1,4): 
                newt.vecs[:,ic] = tet.vecs[:,triad[ic-1],it-i] #three corners for new tet           
            va = newt.vecs[:,1]-newt.vecs[:,0]
            vb = newt.vecs[:,2]-newt.vecs[:,0]
            vc = newt.vecs[:,3]-newt.vecs[:,0]
            newt.vol = abs(dot(va,cross(vb,vc))/2.0)
            newt.fwgt = tet.fwgt[it-i,:]*newt.vol/tet.vol[it-i]
            newt.vwgt =  tet.vwgt[it-i]*newt.vol/tet.vol[it-i]          
            newt.idk = sort(append(array([kp.nk]),[tet.idk[ic,it-i] for ic in triad])) 
            tet.add(newt)
        tet.remove(it-i)#remove old tetrahedron.  it-i below will decrement the values in tetpartial as we delete tets
        kp.make_ids(tet) #easier than trying to change it. 
       # add new tets
        for j in range(tet.ntet):
            if j < 4: #this is connected to only 4
                newk.idt[j] = tet.ntet+j
            else:
                newk.idt[j] = -1
        kp.add(newk)
 
class query:
    def __init__(self):
        ''''''         
    def sum_fermi(self,tet,kp):
        eps = 1e-6 
        sumfwgt = 0.0
        sumpartial = 0.0
        sumfilled = 0.0
        for it in range(tet.ntet):
            for ib in range(kp.nbands):
                sumfwgt += tet.fwgt[it,ib]
                if tet.fwgt[it,ib]>eps and abs(tet.fwgt[it,ib] - tet.fwgt[it,0])>eps:
    #                print it,ib,tet.fwgt[it,ib]
                    sumpartial += tet.fwgt[it,ib]
                else:
                    sumfilled  += tet.fwgt[it,ib]            
    
        print 'fermi weight of partial tets', sumpartial
        print 'fermi weight of filled tets', sumfilled
        print 'total fermi weight', sumfwgt   
        return sumfwgt
    
    def tetwgt(self,tet,kp):
        ''' '''
        print 'Tet, band, weight' 
        for it in range(tet.ntet):
           print [it,[tet.fwgt[it,ib]  for ib in range(kp.nbands)]]

    def ener(self,kp):
        ''' '''
        print 'Kpt, band, ev' 
        for ik in range(kp.nk):
           print [ik,[kp.ener[ik,ib]  for ib in range(kp.nbands)]]
           
#    def ebands(self,tet,kp):
#        ''''''
#        etot = 0.0
#        for it in range(tet.ntet):
#            for ib in range(kp.nbands):
#                for ic in range(4):
#                    etot += tet.fwgt[it,ib]
            
            
#    
#    def tetwgt2kpwgt(self,tet,kp):
#        kw = zeros(kp.nk,dtype=float)
#        for ik in range(kp.nk):
#            for it in range(tet.ntet):
#                for ic in range(4):
#                    if ik == tet.idk[ic,it]: 
#                        kw[ik] += tet.vwgt[it]/4.0
#        print 'Unscaled k point weights from tets'
#        for ik in range(kp.nk): print ik, kw[ik]
#        print 'Total weight',sum(kw)
#                                     
                                                  
################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alal/test/f1_2/'
maindir = '/fslhome/bch/cluster_expansion/alal/test/f1_4/'
tet = tetrahedrons()
kp = kpoints()
wrtvasp = write_vasp()
q = query()

set_printoptions(precision=12)
os.chdir(maindir)
dir = '1init/'
tet.get_tetweights(dir)

kp.get_eigen(dir)
q.ener(kp)
tet.get_tetvecs(dir)
#q.tetwgt(tet,kp) #prints

#print 'tet8', sum(tet.fwgt[8,:]);tet8=sum(tet.fwgt[8,:])
print 'tet vecs'
print tet.vecs
tet.rd_ibzkpt(dir)
print 'tet ids to kpoints'
kp.twgt2kwgt(tet) #kpweights from tetweights
print tet.idk
print 'tet volume weights'
print tet.vwgt
print 'Sum volume weights', sum(tet.vwgt)
print 'tet volumes'
print tet.vol
kp.rd_ibzkpt(dir)
print 'kpoint vecs'
print kp.vecs
print 'nk',kp.nk
kp.make_ids(tet)
print 'kpoint connection to tets'
print kp.idt

sumfwgt1 = q.sum_fermi(tet,kp)
add_tets(tet,kp)

sumfwgt2 = q.sum_fermi(tet,kp)
kp.twgt2kwgt(tet) #kpweights from tetweights
#print 'tet weights'
#print tet.ntet,tet.fwgt
#print 'Energies'
#print kp.ener
#print 'tet vecs'
#print tet.vecs
#print 'tet volume weights'
#print tet.vwgt
#print 'tet volumes'
#print tet.vol
#print 'tet ids to kpoints'
#print tet.idk
#print 'kpoint vecs'
#print kp.vecs
#print 'kpoint to tet ids'
#kp.make_ids(tet)
#print kp.idt
#print 'nk',kp.nk
#print 'ntet',tet.ntet
#print 'Sum volume weights', sum(tet.vwgt)

#### VASP run with added k's
dir = '2addk/'
wrtvasp.kpfile(tet,kp,dir)  
#run vasp !!

print 'Get new eigenvals, weights, but not vecs'
tet.get_tetweights(dir)

kp.twgt2kwgt(tet) #kpweights from tetweights
kp.get_eigen(dir)
q.ener(kp)

sumfwgt3 = q.sum_fermi(tet,kp)


#print 'tet 8old', tet8
#print 'tet 26-29' ,sum(tet.fwgt[26:30,:])
#print sumfwgt3
if abs(sumfwgt3-sumfwgt1)>1e-7: print '  SOMEHOW LOSING FERMI WEIGHT IN VASP!'    

#### Add more kpoints
##weights are getting below integers:
#
add_tets(tet,kp)
tet.volscale = tet.volscale/4.0
tet.vwgt = tet.vwgt*4.0
#
##### VASP run with added k's
dir = '3addk/'
wrtvasp.kpfile(tet,kp,dir)   



print 'Done'

