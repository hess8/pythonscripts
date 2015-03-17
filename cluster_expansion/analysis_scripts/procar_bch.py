#!/usr/bin/python
'''    
Written by Bret Hess, Brigham Young University, 2015
'''
import sys,os,subprocess
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from numpy import zeros,transpose,array,sum,float64,rint,mean,set_printoptions,s_,\
    delete,append,cross,dot,nditer,exp,cumsum
from numpy.linalg import norm
#from numpy.ndarray import flatten
from analysisToolsVasp import nstrip, getEf, readfile

from plotTools import plotxy
from pylab import *
from copy import deepcopy
from itertools import chain

class procar:
    def __init__(self):
#        self.path = []
    
        self.nk = 0
        self.kvecs = []
        self.kweights = []
        self.nbands = 0
        self.weights = []
        self.nspin = 1 #1: no vasp spin, 2: spin
        self.orbs = []
        self.norbs = 0
        self.nions = 0
        self.ener = []
        self.occ = []
        
    def incar(self,dir):
        print 'reading ISPIN and LORBIT from INCAR'
        lines = nstrip(readfile(dir+'INCAR'))
        for line in lines:
            if ('ISPIN' in line) or ('ispin' in line): self.nspin = int(line.split('=')[1])
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
        self.ener = zeros((self.nk, self.nbands, self.nspin),dtype = float) 
        self.occ = zeros((self.nk, self.nbands,self.nspin),dtype = float) 
        self.weights = zeros((self.nk,self.nbands,self.nions,len(self.orbs),self.nspin),dtype = float)
        
        #read weights
        ipos = 2
        for ispin in range(self.nspin):
            ipos += 1; # print 'ipos for spin %i:  %i' % (ispin,ipos)
            for ik in range(self.nk):
                self.kvecs[:,ik] = lines[ipos].split()[3:6]
                self.kweights[ik] = float(lines[ipos].split()[8])
                ipos += 2
                for ib in range(self.nbands):
                    self.ener[ik,ib,ispin] = float(lines[ipos].split()[4])
                    self.occ[ik,ib,ispin] = float(lines[ipos].split()[7])
                    ipos += 3
                    for i in range(self.nions):
                        if mod(ipos,1e4)==0: print '-', #progress bar
                        tempw = lines[ipos].split()[1:self.norbs+1]
                        self.weights[ik,ib,i,:,ispin] = tempw
                        ipos += 1
                    ipos += 2
                ipos += 1
        print       
        
class dos:
    def __init__(self):     
        ''' Creates a density of states plot from procar data '''
        self.title = []
        self.xlabel = '$E-E_F$ (eV)'
        self.ylabel = 'Density of states (arb. units)'
        self.legend = []
        self.colors = []
        self.fig = []
        self.ax1 = []
#        self.maxY =0 #the largest y in the region of interest
        self.minEplot = [] #the lowest energy to include in the plot.  Then maxY will be found in this range...
       
    def plotinit(self,title1,minEplot):
        self.minEplot = minEplot
        self.fig = figure()
        #rcParams['axes.color_cycle']=['r','g']
        self.ax1 = self.fig.add_subplot(111)
        self.ax1.set_color_cycle(['r','b','g', 'k','c', 'm', 'y'])
        title(title1)
        xlabel(self.xlabel)
        ylabel(self.ylabel)  

#        for i in range(N):
#            ax1.loglog(NkfullBZ, deviations[:,i],label=enerlabels[i],linestyle='None',color=cm.jet(1.*(i+1)/N), marker = 'o') # marker = 'o',
    
    def plotline(self,procar,pr_slice,Ef,leg_in):
        '''Adds a line to a plot'''
        [e_arr,d_arr] = self.plotcalc(procar,pr_slice,Ef)                            
        plot(e_arr,d_arr,label=leg_in)

    def plotcumsumline(self,procar,pr_slice,Ef,leg_in):
        '''Adds a line to a plot, just as with plotline, but the only difference is 
        that is plots the cumulutive sum'''
        [e_arr,d_arr] = self.plotcalc(procar,pr_slice,Ef)    
        plot(e_arr,cumsum(d_arr),label=leg_in)

    def plotcumEsumline(self,procar,pr_slice,Ef,leg_in):
        '''Adds a line to a plot, just as with plotline, but the only difference is 
        that is plots the cumulative total E: sum(n(E)E)''' 
        [e_arr,d_arr] = self.plotcalc(procar,pr_slice,Ef)
        temp_arr = e_arr * d_arr  #element wise mult 
        plot(e_arr,cumsum(temp_arr),label=leg_in)

    def plotcalc(self,procar,pr_slice,Ef):
        dE = 0.03 #plot resolution (bin width in eV
        ngauss = 10 #gaussian width is 2ngauss+1 plot bins
        eshifted = procar.ener - Ef
        emin = amin(eshifted.flatten()) 
        emax = amax(eshifted.flatten())
        imin = max(0,self.iener(self.minEplot,emin,ngauss,dE)) 
        nE = int(ceil((emax - emin)/dE)) + 1 + 2*ngauss
        d_arr = zeros(nE,dtype = float)
        e_arr = zeros(nE,dtype = float)
        e_arr = [emin-ngauss*dE + i*dE for i in range(nE)]
        gaussw = array([exp(-(float(i)/float(ngauss/2))**2) for i in range(-ngauss,ngauss+1)])       
        gaussw = gaussw/sum(gaussw) #normalized to 1
        sum1 = 0.0
        for ispin in self.sliceparse(pr_slice.spin,procar.nspin):
            for ik in self.sliceparse(pr_slice.ks,procar.nk):
                kw = procar.kweights[ik]
                for ib in self.sliceparse(pr_slice.bands,procar.nbands):
                    wener = 0.0
                    occ = procar.occ[ik,ib,ispin]
                    for ii in self.sliceparse(pr_slice.ions, procar.nions):
                        for io in self.sliceparse(pr_slice.orbs,procar.norbs):                     
                            wener += procar.weights[ik,ib,ii,io,ispin]   
                            sum1 += procar.weights[ik,ib,ii,io,ispin]        
                    wener = wener * kw #* occ  #leave off the factor occ if you want to plot empty states!
                    if wener > 0.0:
                        ie = self.iener(eshifted[ik,ib,ispin],emin,ngauss,dE)
                        d_arr[ie-ngauss:ie+ngauss+1] += wener*gaussw 
        print 'integrated electron weight in line', sum1  
        return e_arr[imin:],d_arr[imin:]              
       
    def iener(self,e,emin,ngauss,dE):
        '''Gives closest index in e_array for a particular energy.  emin maps has position ngauss above zero'''
        return int(ceil((e-emin)/dE ))+ ngauss
        
    def sliceparse(self,slicelist,N):
        '''Slices here are just lists of indices. Must convert from VASP's numbering to python's (-1)
        2:4 includes 2,3,4'''
        list1 = []
        for slicetxt in slicelist:
            if '-' in slicetxt:
                list1.append(range(int(slicetxt.split('-')[0]),int(slicetxt.split('-')[1])))
            elif ':' in slicetxt: # e.g.  "6:" input is "5:" output  ":7" input is ":7" output'
                splt = slicetxt.split(':')
                if splt[0] == ':' : list1.append(range(0,int(splt[1])))
                else: list1.append(range(int(splt[0])-1,N))
            elif slicetxt == 'all':
                list1.append(range(N))
            else: 
                list1.append([int(slicetxt)-1])
        return list(chain(*list1)) #flattens list
        
 #        '''slice function is in form: slice(start, stop, step) for start:stop:step'''
#        if '-' in slicetxt:
#            return slice(int(slicetxt.split('-')[0]),int(slicetxt.split('-')[0]))
#        elif ':' in slicetxt:
#            splt = slicetxt.split(':').replace('',None)
#            return slice(splt[0],splt[1])
#        elif slicetxt == 'all':
#            return slice(None,None)
#        else: 
#            return slice(int(slicetxt), int(slicetxt)+1)   
    def plotend(self,file):    
        legend(loc='upper left')    
        show()
        self.fig.savefig(file)
        close 
class pr_slice():
    def __init__(self): 
        '''defines what slices of procar weights to add to a single plot line''' 
        self.ions = ['all']
        self.bands = ['all']
        self.ks = ['all']
        self.spin = ['all']
        self.orbs = ['all']
#        dosarr = zeros()                                                                 
################# script #######################

'''Final all atomic folders (capitalized names) in maindir, and identify which stucture folders 
have DOS runs.  Create '''
#maindir = '/fslhome/bch/cluster_expansion/alal/test/f1_2/'
#maindir = '/fslhome/bch/cluster_expansion/graphene/analysis/1220/DOS/'
#maindir = '/fslhome/bch/cluster_expansion/graphene/analysis/1220/DOSlorbit11x30/'
atom1 = 'V'
atom2 = 'C'
maindir = '/fslhome/bch/cluster_expansion/graphene/analysis/{}_sv_simpleDOS/3/DOS/'.format(atom1)

#maindir = '/fslhome/bch/vasprun/graphene.structures/transmet.half_graphane/half_graphane/dos/adatom_W/ISIF_4/IBRION_2/ISPIN_2/dos/'
os.chdir(maindir)
Ef = float(getEf(maindir))
pcr = procar();dos = dos();pr_slice = pr_slice()
pcr.incar(maindir)
pcr.read(maindir+'PROCAR')

#DOS
#dos.plotinit('W Structure 1220')
#pr_slice.ions = ['1-12']
#dos.plotline(pcr,pr_slice,Ef,'C atoms')
#pr_slice.ions = ['13-16']
#dos.plotline(pcr,pr_slice,Ef,'H atoms')
#pr_slice.ions = ['17:']
#dos.plotline(pcr,pr_slice,Ef,'W atoms')
#pr_slice.ions = ['all']
#dos.plotline(pcr,pr_slice,Ef,'Total')

minEplot = -100 #lower x axis limit
dos.plotinit('DOS {} Structure 3'.format(atom1),minEplot) #title
pr_slice.ions = ['3:'] #atoms start counting at 1
pr_slice.orbs = ['1']
dos.plotline(pcr,pr_slice,Ef,'S of {}'.format(atom1))
pr_slice.orbs = ['2-4']
dos.plotline(pcr,pr_slice,Ef,'P of {}'.format(atom1))
pr_slice.orbs = ['5-9']
dos.plotline(pcr,pr_slice,Ef,'D of {}'.format(atom1))
pr_slice.orbs = ['all']
dos.plotline(pcr,pr_slice,Ef,'Total {}'.format(atom1))

pr_slice.ions = ['1-2'] #atoms start counting at 1
pr_slice.orbs = ['1']
dos.plotline(pcr,pr_slice,Ef,'S of C')
pr_slice.orbs = ['2-4']
dos.plotline(pcr,pr_slice,Ef,'P of C')
pr_slice.orbs = ['all']
dos.plotline(pcr,pr_slice,Ef,'Total C')

dos.plotend('DOS_{}C_orbs'.format(atom1)) #file 
#=======================
# Make cumulative sum plots
dos.plotinit('Cumulative DOS {} Structure 3'.format(atom1),minEplot) #title
pr_slice.ions = ['3:'] #atoms start counting at 1
pr_slice.orbs = ['1']
dos.plotcumsumline(pcr,pr_slice,Ef,'S of {}'.format(atom1))
pr_slice.orbs = ['2-4']
dos.plotcumsumline(pcr,pr_slice,Ef,'P of {}'.format(atom1))
pr_slice.orbs = ['5-9']
dos.plotcumsumline(pcr,pr_slice,Ef,'D of {}'.format(atom1))
pr_slice.orbs = ['all']
dos.plotcumsumline(pcr,pr_slice,Ef,'Total {}'.format(atom1))

pr_slice.ions = ['1-2'] #atoms start counting at 1
pr_slice.orbs = ['1']
dos.plotcumsumline(pcr,pr_slice,Ef,'S of C')
pr_slice.orbs = ['2-4']
dos.plotcumsumline(pcr,pr_slice,Ef,'P of C')
pr_slice.orbs = ['all']
dos.plotcumsumline(pcr,pr_slice,Ef,'Total C')

dos.plotend('CumDOS_{}C_orbs'.format(atom1)) #file 
       
# Make cumulative total energy plots
dos.plotinit('Total energy {} Structure 3'.format(atom1),minEplot) #title
pr_slice.ions = ['3:'] #atoms start counting at 1
pr_slice.orbs = ['1']
dos.plotcumEsumline(pcr,pr_slice,Ef,'S of {}'.format(atom1))
pr_slice.orbs = ['2-4']
dos.plotcumEsumline(pcr,pr_slice,Ef,'P of {}'.format(atom1))
pr_slice.orbs = ['5-9']
dos.plotcumEsumline(pcr,pr_slice,Ef,'D of {}'.format(atom1))
pr_slice.orbs = ['all']
dos.plotcumEsumline(pcr,pr_slice,Ef,'Total {}'.format(atom1))

pr_slice.ions = ['1-2'] #atoms start counting at 1
pr_slice.orbs = ['1']
dos.plotcumEsumline(pcr,pr_slice,Ef,'S of C')
pr_slice.orbs = ['2-4']
dos.plotcumEsumline(pcr,pr_slice,Ef,'P of C')
pr_slice.orbs = ['all']
dos.plotcumEsumline(pcr,pr_slice,Ef,'Total C')

dos.plotend('CumETot_{}C_orbs'.format(atom1)) #file 
print 'Done'
