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
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
sys.path.append('/fslhome/bch/graphener/graphener')
from kmeshroutines import svmesh,  svmeshNoCheck,svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection


class meshConstruct(): 
    '''Makes superperiodic structures and cuts'''
    from comMethods import readfile,writefile,trimSmall,areEqual
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform

    def __init__(self):
        '''init'''
         

    def relaxMeshSym(self,reciplatt,targetNmesh,path):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        nCoarse = targetNmesh #for now
        self.subRandom(nCoarse,reciplatt)
        for i in range(nCoarse):
            print self.pos[i]
        self.plot2dPts('x','y') 
        self.plot2dPts('x','z')          
        sys.exit('stop')
        
        
        return meshvecs, Nmesh, lattype, pfB, pf, status
    
    def subRandom(self,nCoarse,B):
        self.pos = zeros((nCoarse,3))
        self.subrand = rand(3)
#         ns = svmeshNoCheck(nCoarseB)
        afact = sqrt(2)
        bfact = sqrt(3)
        cfact = sqrt(5)
        for ipos in range(nCoarse):
            self.subrand = mod(self.subrand + array([afact,bfact,cfact]), array([1,1,1]))
            self.pos[ipos,:] = dot(B,self.subrand)
    
    def plot2dPts(self,ax0,ax1):
        fig = figure()
        i0 = ('x','y','z').index(ax0)
        i1 = ('x','y','z').index(ax1)
        scatter(self.pos[:,i0],self.pos[:,i1])
        xlabel('{}'.format(ax0))
        ylabel('{}'.format(ax1))
        name = '{}-{}'.format(ax0,ax1)
        title(name)
        fig.savefig(name+'.pdf');
        show()
            

                    
                
            
            