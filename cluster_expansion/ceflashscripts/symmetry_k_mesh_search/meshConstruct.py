'''
Builds a mesh from the symmetry operations and random starting points.
There is no lattice
0. The mesh should be optimized as an atomic space superlattice (transform of k-space)
0.1 Choose a 
1. nCoarse sub-random points in the BZ cell parallelpiped.  
    We don't need to work in Voronoi cell version because we are running in 3x3x3 supercell of the BZ
    https://en.wikipedia.org/wiki/Low-discrepancy_sequence#Additive_recurrence for subrandom 
2. Define cutoff as rc = (nCoarse/Vcell)^1/3
3. Make all copies due to translation and symmetry operations
4. In neighboring cells (3x3x3 construction of plpd), keep only those that 
are within e.g. 1.5 * rc of the cell boundaries
5. If two points are less than e.g. 0.8*rc away from each other, combine them 
into a point at their center. 
5.5 if we have too many delete those that are closest to their neighbors
5.6 If we have too few, get the Delaunay tetrahedra, and add to the center 
of the largest tets. 
6. Keep track of independent points and symmetry partners
Ind. points are in the original unit cell, closest to the origin, smallest x comp, smallest y comp, smallest z
to break ties.
7. Calculate the force on each ind. point from all the points that lie within
a distance rForce.  Use a power law: Fi = k*Sum(r_vec_i,j/ rij^p
8.  Move the ind. points in a vector proportional to the force
9.  Copy all the ind. points by trans and pt symmetry again. 
10. repeat from 7 until the force on each point is less than minForce
11. Copy all the ind. points by trans and pt symmetry again. 
________
Fine tetrahedral mesh
________
1. Find Delaunay tetrahedra
2. Subdivide each tet that has an independent point in it with a formula that 
minimizes the S/V of the 6 new tets. Call these independent children, 
even though they might have parents that are not independent
3. Copy these by trans and point sym, removing duplicates
4. Subdivide the indChildren (repeat 2 and 3) until we have the number of points we want.
5. By symmetry, find the number of truly independent among all points, and their weights. 
6. Write these to KPOINTS: https://cms.mpi.univie.ac.at/vasp/vasp/Entering_all_k_points_explicitly.html 
until we have the number of points we want. 

All matrices store vectors as COULUMNS
'''
import os, subprocess,sys,re,time
from numpy import mod,dot, transpose, rint,floor,ceil,zeros,array,sqrt,average,std,amax,amin,int32,sort,count_nonzero,\
    delete,mean,square,argmax,argmin,insert
from numpy.linalg import inv, norm, det
from numpy.random import rand
from copy import copy,deepcopy
from sched import scheduler
from itertools import chain
from matplotlib.pyplot import subplots,savefig,imshow,close,plot,title,xlabel,ylabel,figure,show,scatter
import matplotlib.image as mpimg
import datetime
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
sys.path.append('/fslhome/bch/graphener/graphener')
from kmeshroutines import svmesh,  svmeshNoCheck,svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection

def timestamp():
    ts = time.time()
    return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

class meshConstruct(): 
    '''Makes superperiodic structures and cuts'''
    from comMethods import readfile,writefile,trimSmall,areEqual
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform

    def __init__(self):
        '''init'''
         

    def relaxMeshSym(self,B,targetNmesh,path):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        self.nCoarse = targetNmesh
        vol = det(B)
        self.edgeFactor = 1.5
        rc = self.edgeFactor*(targetNmesh/vol)^(1/3.0) #cutoff for forces, neighbors in other cells. 
        self.initPoints(B) 
        [self.symops,self.nops] = getGroup(B)
        print 'Number of symmetry operations:', self.nops
        self.transPoints
        self.symPoints()
               
        return meshvecs, Nmesh, lattype, pfB, pf, status
    
    def transPoints():
        '''Make copies of points in 3x3x3 supercell, but keep only those that 
        are within rc of edges'''
#         self.posAllCells = zeros((1,1,1,self.nCoarse,3)) #first three digits will be a,b,c multiples to get the 27 cells:  -1, 0, 1
        self.nextLabel = 0
        self.points = zeros((len(self.nCoarse*27),dtype = [('label', 'int8'),('dep', 's1'),('sponsor', 'int8'),('pos', '3float'),('cell', '3int')])
        aEdge = rc/norm(B[:,0]); bEdge = rc/norm(B[:,1]); cEdge = rc/norm(B[:,2])
#         bigDirCell = array([-1-])
        ipoint = 0
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    for ipos in range(self.nCoarse):
                        dpos = self.directPos[ipos]
                        pos = self.pos[ipos]
                        newDirPos = dpos + array([i,j,k])
                        if  (-1-aEdge < dpos[0] < 1 + aEdge) and (-1-bEdge < dpos[2] < 1 + bEdge) and (-1-cEdge < dpos[2] < 1 + cEdge):
                            newPos = pos + i*B[:,0] + j*B[:,1] + k*B[:,2]
                            self.points[ipoint]['label'] = self.nextLabel
                            self.points[ipoint]['dep'] = 'U' #unknown
                            self.points[ipoint]['sponsor'] = -1 #not assigned
                            self.points[ipoint]['cell'] = [i,j,k]
                            self.points[ipoint]['pos'] = newPos
                            ipoint += 1
        self.nextLabel = ipoint     
    
    def initPoints(self,B):
        self.subRandom(B)  
        self.plotPos()      
        sys.exit('stop')
            
    def subRandom(self,B):
        self.pos = zeros((self.nCoarse,3))
        self.directPos = zeros((self.nCoarse,3))
        self.subrand = rand(3)
        afact = sqrt(2)
        bfact = sqrt(3)
        cfact = sqrt(5)
        for ipos in range(self.nCoarse):
            self.subrand = mod(self.subrand + array([afact,bfact,cfact]), array([1,1,1]))
            self.directPos[ipos,:] = self.subrand
#             print 'srand',self.subrand
            self.pos[ipos,:] = dot(transpose(B),self.subrand)
    
    def plotPos(self):
        self.plot2dPts(timestamp,'x','y') 
        self.plot2dPts(timestamp,'x','z')  
        
    def plot2dPts(self,tag,ax0,ax1):
        fig = figure()
        i0 = ('x','y','z').index(ax0)
        i1 = ('x','y','z').index(ax1)
        scatter(self.pos[:,i0],self.pos[:,i1])
        xlabel('{}'.format(ax0))
        ylabel('{}'.format(ax1))
        name = '{} {}-{}'.format(tag,ax0,ax1)
        title(name)
        fig.savefig(name+'.pdf');
        show()
