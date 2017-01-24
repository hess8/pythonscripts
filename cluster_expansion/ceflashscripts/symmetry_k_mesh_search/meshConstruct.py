'''
Builds a mesh from the symmetry operations and random starting points.
There is no lattice
0. The mesh should be optimized as an atomic space superlattice (transform of k-space)
0.1 Choose a 
1. nCoarse sub-random points in the BZ cell parallelpiped.  
    We don't need to work in Voronoi cell version because we are running in 3x3x3 supercell of the BZ
    https://en.wikipedia.org/wiki/Low-discrepancy_sequence#Additive_recurrence for subRandSym 
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
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection,intoVoronoi,intoCell

def timestamp():
    return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())

class meshConstruct(): 
    '''Makes superperiodic structures and cuts'''
    from comMethods import readfile,writefile,trimSmall,areEqual,directFromCart,cartFromDirect
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform

    def __init__(self):
        '''init'''
         
    def relaxMeshSym(self,B,targetNmesh,path):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        self.B = B
        [self.symops,self.nops] = getGroup(self.B)
        self.nCoarse = targetNmesh
        self.path = path
        print 'Number of desired points:', targetNmesh
        print 'Symmetry operations:', self.nops
        vol = abs(det(B))
        ravg = (vol/targetNmesh)**(1/3.0)
        self.edgeFactor = 1.5
        self.rmin = 0.8*ravg #
        self.rmax = self.edgeFactor*ravg #cutoff for forces, neighbors in other cells. 
        print 'Cutoff rc:', self.rmax
        self.initPoints() 
        self.transPoints()
        self.relax()
        
        sys.exit('stop')     
        return meshvecs, Nmesh, lattype, pfB, pf, status
       
    def relax(self):
        '''Conjugate gradient relaxation of the potential energy'''
        
    def energy(self):
        '''restrict the energy sum to the pairs that contain independent points'''
        ener = 0.0
        p = 4
        for ipos in self.ind:
            for jpos in range(self.nCoarse):
                if ipos != jpos:
                    r = norm(self.points[ipos]['pos']-self.points[jpos]['pos'])
                    ener += 1/r**p    

    def transPoints(self):
        '''Make copies of points in 3x3x3 supercell, but keep only those that 
        are within rc of edges'''
#         self.posAllCells = zeros((1,1,1,self.nCoarse,3)) #first three digits will be a,b,c multiples to get the 27 cells:  -1, 0, 1
#         self.points = zeros((len(self.nCoarse*27),dtype = [('label', 'int8'),('dep', 's1'),('sponsor', 'int8'),('pos', '3float'),('cell', '3int')])
        self.points = zeros(self.nCoarse*27,dtype = [('label', 'int8'),('dep', 'S1'),('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('cell', '3int')])
#         self.labels = [] # list of active points labels 

        aEdge = self.rmax/norm(self.B[:,0]); bEdge = self.rmax/norm(self.B[:,1]); cEdge = self.rmax/norm(self.B[:,2])
#         bigDirCell = array([-1-])
        
#         self.points[:self.nCoarse] = self.cellPoints[:self.nCoarse]
        ipoint = 0
        for i in [0,-1,1]:
            for j in [0,-1,1]:
                for k in [0,-1,1]:
                    for ipos in range(self.nCoarse):
                        dpos = self.cellPoints[ipos]['dpos'] #dpos will always be the direct coordinates in the original cell
                        #move to parallepiped cell
                        pos = dot(transpose(self.B),dpos)
                        self.cellPoints[ipos]['pos']  = pos #was voronoi position before                
                        newDirPos = dpos + array([i,j,k])
                        if  (-1-aEdge < newDirPos[0] < 1 + aEdge) and (-1-bEdge < newDirPos[2] < 1 + bEdge) and (-1-cEdge < newDirPos[2] < 1 + cEdge):
                            newPos = dot(transpose(self.B),newDirPos)
#                             self.points[ipoint]['label'] = ipoint
                            self.labels.append(ipoint)
                            if i != 0 or k !=0 or j!=0: #then this point is dependent by translation
                                self.points[ipoint]['dep'] = 'T'
                                self.points[ipoint]['sponsor'] = ipos #not assigned
                                self.dep.append(ipoint)                        
                            self.points[ipoint]['cell'] = [i,j,k]
                            self.points[ipoint]['pos'] = newPos
                            self.points[ipoint]['dpos'] = newDirPos
                            ipoint += 1
        self.npoints = ipoint-1
        self.plotPos(self.points,'pos')
        self.plotPos(self.points,'dpos')
        return   
    
    def initPoints(self):
        self.subRandSym()  
#         self.combineNear() #caution, this does not preserve symmetry
#         self.plotPos(self.cellPoints)  
        
#     def combineNear(self): 
#         '''Doesn't preserve symmetry, so not using it yet'''
#         self.cellPoints2 = zeros(self.nCoarse+self.nops,dtype = [('label', 'int8'),('dep', 'S1'),('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('cell', '3int')])
# #         combine = []
#         for jpos in range(self.nCoarse):
#             for kpos in range(jpos+1,self.nCoarse):
#                 r = norm(self.cellPoints[jpos]['pos'] - self.cellPoints[kpos]['pos']) < self.rmin
#                 if 0.0 < norm(self.cellPoints[jpos]['pos'] - self.cellPoints[kpos]['pos']) < self.rmin:
#                     print jpos,kpos
#                     print 'r', norm(self.cellPoints[jpos]['pos'] - self.cellPoints[kpos]['pos'])
#                     if kpos in self.labels: self.labels.remove(kpos)
#                     if kpos in self.ind: self.ind.remove(kpos)
#                     if kpos in self.dep: self.dep.remove(kpos)
#                     self.cellPoints[jpos]['pos'] = (self.cellPoints[jpos]['pos'] + self.cellPoints[kpos]['pos'])/2.0     
#         ipos = 0
#         for label in self.labels:
#             self.cellPoints2[ipos] = self.cellPoints[label]
#             ipos += 1
#         self.nCoarse = ipos
#         self.cellPoints[label] = self.cellPoints2[ipos]
#         print 'Points generated:',self.nCoarse
#         return
                      
    def subRandSym(self):
        self.labels = [] # list of active points labels
        self.cellPoints = zeros(self.nCoarse+self.nops,dtype = [('label', 'int8'),('dep', 'S1'),('sponsor', 'int8'),('pos', '3float'),('vpos', '3float'),('dpos', '3float'),('cell', '3int')])
        self.directPos = zeros((self.nCoarse,3))
        self.subrand = rand(3)
        afact = sqrt(2)
        bfact = sqrt(3)
        cfact = sqrt(5)
        self.ind = []
        self.dep = []
        ipos = -1
#         for i in range(self.nops):
#             print i
#             print self.symops[:,:,i]
        while ipos < self.nCoarse:
            ipos += 1
            self.subrand = mod(self.subrand + array([afact,bfact,cfact]), array([1,1,1]))
            temp = self.subrand
            pos = dot(transpose(self.B),temp)
            pos = intoVoronoi(pos,self.B)#now in Voronoi cell
            self.cellPoints[ipos]['dep'] = 'I' #unknown
            self.cellPoints[ipos]['sponsor'] = ipos                              
            self.cellPoints[ipos]['pos'] = pos 
            self.cellPoints[ipos]['dpos'] = temp #dpos will always be the direct coordinates in the original cell
            self.labels.append(ipos)
            self.ind.append(ipos)
            iInd = ipos
            #now all symmetry operations will keep the point in the BZ 
            sympoints = [str(pos[0])[:5]+str(pos[1])[:5]+str(pos[2])[:5]]
            for op in range(self.nops):
                newpos = dot(self.symops[:,:,op],pos)
                posStr = str(newpos[0])[:7]+str(newpos[1])[:7]+str(newpos[2])[:7]
                if posStr not in sympoints:
                    sympoints.append(posStr)
                    ipos += 1
                    self.labels.append(ipos)
                    self.dep.append(ipos)
                    self.cellPoints[ipos]['dep'] = 'S'
                    self.cellPoints[ipos]['sponsor'] = iInd 
#                     newpos=intoVoronoi(newpos,self.B)                                
                    self.cellPoints[ipos]['dpos'] = self.directFromCart(self.B,intoCell(newpos,self.B)) 
                    checkpos = intoVoronoi(newpos,self.B)
                    if norm(checkpos-newpos)> 1e-6:
                        print 'New Voronoi cell pos for',newpos,checkpos                   
                    self.cellPoints[ipos]['pos'] = newpos    
        self.nCoarse = ipos  
        self.plotPos(self.cellPoints,'pos')
#         self.plotPos(self.cellPoints,'dpos')
        print 'Points in unit cell:',self.nCoarse                           
        return
    
#     def tooNear(pos,ipos):
#         for testpos in self.points[:ipoint]['pos']:
#             if norm(pos-testpos) < self.rmin:
#                 return False
#         return True
    
    def plotPos(self,arr,field):
        self.plot2dPts(arr,field,timestamp(),'x','y') 
        self.plot2dPts(arr,field,timestamp(),'x','z')  
        
    def plot2dPts(self,arr,field,tag,ax0,ax1):
        fig = figure()
        i0 = ('x','y','z').index(ax0)
        i1 = ('x','y','z').index(ax1)
        scatter(arr[:self.nCoarse][field][:,i0],arr[:self.nCoarse][field][:,i1])
        xlabel('{}'.format(ax0))
        ylabel('{}'.format(ax1))
        name = '{}_{}_{}-{}'.format(tag,field,ax0,ax1)
        title(name)
        fig.savefig(self.path+'/'+ name+'.pdf');
        os.system('cp {} {}'.format(self.path+'/'+ name+'.pdf',self.path+ '/latest{}-{}'.format(ax0,ax1)+'.pdf' ))
        show()
