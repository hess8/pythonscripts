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
6. Keep track of independent points (all are "sponsors") and symmetry partners
Dependent points (which can also be sponsors) have "dep" label:  if by translations: "-1,0,1". 
if by symmetry, a single positive integer "3" that gives the operation
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
# from scipy.optimize import fmin_cg
from numpy.linalg import inv, norm, det
from numpy.random import rand
from copy import copy,deepcopy
from sched import scheduler
from itertools import chain
from matplotlib.pyplot import subplots,savefig,imshow,close,plot,title,xlabel,ylabel,figure,show,scatter
import matplotlib.image as mpimg
import datetime
from _ast import operator
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
sys.path.append('/fslhome/bch/graphener/graphener')


# from conjGradMin.optimize import fmin_cg

from kmeshroutines import svmesh,  svmeshNoCheck,svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection,intoVoronoi,intoCell

def timestamp():
    return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())

# def energy(indComps,parent):
#     '''restrict the energy sum to the pairs that contain independent points'''
#     indVecs = parent.indComps.reshape((len(parent.ind),3))
#     #update all positions by symmetry and translation
#     parent.oldPoints = parent.points
#     parent.updatePoints(indVecs)
# #         parent.plotPos(parent.points,parent.npoints,'pos')
# #         for i in range(20):
# #             print 'old',parent.oldPoints[i]['pos']
# #             print 'new',parent.points[i]['pos']
# #             print
#     for i in range(parent.npoints):
#         move = norm(parent.points[i]['pos']-parent.oldPoints[i]['pos'])
#         if move>0:
#             print i,move
#     ener = 0.0
#     p = 4.0
#     scale = 1
#     for ipos in range(len(parent.ind)):
#         for jpos in range(parent.npoints):
#             r = norm(indVecs[ipos]-parent.points[jpos]['pos'])
#             if r > 1e-4*parent.rmin:
#                 ener += scale*(1/r)**p
# #                 
#     print 'energy:',parent.count, ener
# #     if parent.count == 19:
# #         parent.plotPos(parent.points,parent.npoints,'pos')
#         
#     parent.count += 1
#     #now update all the dependent positions
#   #  ! Note:  need to update the positions at each step!  Do we have to get inside
#     return ener


class meshConstruct(): 
    '''Makes superperiodic structures and cuts'''
    from comMethods import readfile,writefile,trimSmall,areEqual,directFromCart,cartFromDirect
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
    from conjGradMin2 import fmin_cg,minimize_cg,approx_fprime

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
        self.ravg = (vol/targetNmesh)**(1/3.0)
        self.edgeFactor = 1.5
        self.rmin = 0.8*self.ravg #
        self.rmax = self.edgeFactor*self.ravg #cutoff for forces, neighbors in other cells. 
        print 'Cutoff rc:', self.rmax
        self.initPoints() 
        self.transPoints()
        self.relax()
        
        sys.exit('stop')     
        return meshvecs, Nmesh, lattype, pfB, pf, status
       
    def relax(self):
        '''Conjugate gradient relaxation of the potential energy.
        fmin_cg(f, x0, options...)
        f : callable, f(x, *args) Objective function to be minimized. Here x must be a 1-D array of the variables that are to be changed in the search for a minimum, and args are the other (fixed) parameters of f.
        x0 : ndarray A user-supplied initial estimate of xopt, the optimal value of x. It must be a 1-D array of values.
 
         >>> args = (2, 3, 7, 8, 9, 10)  # parameter values
        >>> def f(x, *args):
        ...     u, v = x
        ...     a, b, c, d, e, f = args
        ...     return a*u**2 + b*u*v + c*v**2 + d*u + e*v + f
        >>> def gradf(x, *args):
        ...     u, v = x
        ...     a, b, c, d, e, f = args
        ...     gu = 2*a*u + b*v + d     # u-component of the gradient
        ...     gv = b*u + 2*c*v + e     # v-component of the gradient
        ...     return np.asarray((gu, gv))
        >>> x0 = np.asarray((0, 0))  # Initial guess.
        >>> from scipy import optimize
        >>> res1 = optimize.fmin_cg(f, x0, fprime=gradf, args=args)
 
        The energy must be a function of 1-D inputs, so it will be '''
        self.indComps = array([self.points[i]['pos'] for i in self.ind]).flatten()
        initialGuess = [i+ 0.5*self.rmin*rand() for i in self.indComps]  #just add noise to initial positions
        self.count = 0
        
        epsilon = self.ravg/100
        out = self.fmin_cg(initialGuess,epsilon)
        print out
#         for i in range(27):
#             print i ,self.indComps[i]
#         vectors = self.indComps.reshape((self.npoints,3))
#         for i in range(9):
#             print i ,vectors[i,:]        
        return
    
    def updatePoints(self,indVecs):
        for ipoint in range(self.npoints):
            if ipoint in self.ind:
                place = self.ind.index(ipoint)
                self.points[ipoint]['pos'] = indVecs[place,:]
            else:
                dep = self.points[ipoint]['dep']
                sponsor = self.points[ipoint]['sponsor']
#                 print 'test',ipoint
#                 print 'dep' , dep
#                 print 'spons',sponsor
                if dep == 'I':
                    sys.exit('Stop. Point {} is labeled as independent, but is not in self.ind'.format(ipoint))
                elif ',' not in dep: # symmetry dependent
                    self.points[ipoint]['pos'] = dot(self.symops[:,:,int(dep)],self.points[sponsor]['pos'])
                else: #translation
                    ms = array([int(i) for i in dep.split(',')])
                    svec = self.points[sponsor]['pos']
                    self.points[ipoint]['pos'] = svec + dot(self.B,ms)
 
    def energy(self,indComps):
        '''restrict the energy sum to the pairs that contain independent points'''
        indVecs = indComps.reshape((len(self.ind),3))
        #update all positions by symmetry and translation
        self.oldPoints = self.points
        self.updatePoints(indVecs)
#         self.plotPos(self.points,self.npoints,'pos')
#         for i in range(20):
#             print 'old',self.oldPoints[i]['pos']
#             print 'new',self.points[i]['pos']
#             print
        for i in range(self.npoints):
            move = norm(self.points[i]['pos']-self.oldPoints[i]['pos'])
            if move>0:
                print i,move
        ener = 0.0
        p = 4.0
        scale = 1
        for ipos in range(len(self.ind)):
            for jpos in range(self.npoints):
                r = norm(indVecs[ipos]-self.points[jpos]['pos'])
                if r > 1e-4*self.rmin:
                    ener += scale*(1/r)**p
#                 
        print 'energy:',self.count, ener
        if self.count == 19:
            self.plotPos(self.points,self.npoints,'pos')
             
        self.count += 1
        #now update all the dependent positions
      #  ! Note:  need to update the positions at each step!  Do we have to get inside
        return ener

    def transPoints(self):
        '''Make copies of points in 3x3x3 supercell, but keep only those that 
        are within rc of edges'''
        self.points = zeros(self.nCoarse*27,dtype = [('label', 'int8'),('dep', 'S8'),('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('cell', '3int')])
        aEdge = self.rmax/norm(self.B[:,0]); bEdge = self.rmax/norm(self.B[:,1]); cEdge = self.rmax/norm(self.B[:,2])
        ipoint = 0
        for i in [0,-1,1]:
            for j in [0,-1,1]:
                for k in [0,-1,1]:
                    for ipos in range(self.nCoarse):
                        dpos = self.cellPoints[ipos]['dpos'] #dpos will always be the direct coordinates in the original cell
                        #move to parallepiped cell
                        pos = dot(self.B,dpos)
                        self.cellPoints[ipos]['pos']  = pos #was voronoi position before                
                        newDirPos = dpos + array([i,j,k])
                        if  (-aEdge < newDirPos[0] < 1 + aEdge) and (-bEdge < newDirPos[1] < 1 + bEdge) and (-cEdge < newDirPos[2] < 1 + cEdge):
                            newPos = dot(self.B,newDirPos)
#                             self.points[ipoint]['label'] = ipoint
                            self.labels.append(ipoint)
                            if i == 0 and j == 0 and k == 0: #original cell
                                self.points[ipoint]['dep'] = self.cellPoints[ipoint]['dep'] 
                                self.points[ipoint]['sponsor'] = self.cellPoints[ipoint]['sponsor'] 
                            else: #then this point is dependent by translation
                                self.points[ipoint]['dep'] = '{},{},{}'.format(i,j,k)
                                self.points[ipoint]['sponsor'] = ipos #The sponsor may itself be dependent.  
                                self.dep.append(ipoint)                        
                            self.points[ipoint]['cell'] = [i,j,k]
                            self.points[ipoint]['pos'] = newPos
                            self.points[ipoint]['dpos'] = newDirPos
                            ipoint += 1
        self.npoints = ipoint-1
        self.plotPos(self.points,self.npoints,'pos')
#         self.plotPos(self.points,self.npoints,'dpos')
        return   
    
    def initPoints(self):
        self.subRandSym()  

        
#     def combineNear(self): 
#         '''Doesn't preserve symmetry, so not using it yet'''
#         self.cellPoints2 = zeros(self.nCoarse+self.nops,dtype = [('label', 'int8'),('dep', 'S8'),('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('cell', '3int')])
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
        self.cellPoints = zeros(self.nCoarse+self.nops,dtype = [('label', 'int8'),('dep', 'S8'),\
            ('sponsor', 'int8'),('pos', '3float'),('vpos', '3float'),('dpos', '3float'),('cell', '3int')])
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
            pos = dot(self.B,temp)
            pos = intoVoronoi(pos,self.B)#now in Voronoi cell
            self.cellPoints[ipos]['dep'] = 'I' #independent
            self.cellPoints[ipos]['sponsor'] = ipos                              
            self.cellPoints[ipos]['pos'] = pos 
            self.cellPoints[ipos]['dpos'] = temp #dpos will always be the direct coordinates in the original cell
            self.labels.append(ipos)
            self.ind.append(ipos)
            iInd = ipos
            #now all symmetry operations will keep the point in the BZ 
            strLen = 7
            sympoints = [str(pos[0])[:strLen]+str(pos[1])[:strLen]+str(pos[2])[:strLen]]
            for op in range(self.nops):
                newpos = dot(self.symops[:,:,op],pos)
                posStr = str(newpos[0])[:strLen]+str(newpos[1])[:strLen]+str(newpos[2])[:strLen]
                if posStr not in sympoints:
                    sympoints.append(posStr)
                    ipos += 1
                    self.labels.append(ipos)
                    self.dep.append(ipos)
                    self.cellPoints[ipos]['dep'] = '{}'.format(str(op))
                    self.cellPoints[ipos]['sponsor'] = iInd 
#                     newpos=intoVoronoi(newpos,self.B)                                
                    self.cellPoints[ipos]['dpos'] = self.directFromCart(self.B,intoCell(newpos,self.B)) 
#                     checkpos = intoVoronoi(newpos,self.B)
#                     if norm(checkpos-newpos)> 1e-6:
#                         print 'New Voronoi cell pos for',newpos,checkpos                   
                    self.cellPoints[ipos]['pos'] = newpos    
        self.nCoarse = ipos  
#         self.plotPos(self.cellPoints,self.ncoarse,'pos')
#         self.plotPos(self.cellPoints,self.ncoarse,'dpos')
        print 'Points in unit cell:',self.nCoarse                           
        return
    
#     def tooNear(pos,ipos):
#         for testpos in self.points[:ipoint]['pos']:
#             if norm(pos-testpos) < self.rmin:
#                 return False
#         return True
    
    def plotPos(self,arr,npoints,field):
        self.plot2dPts(arr,npoints,field,timestamp(),'x','y') 
        self.plot2dPts(arr,npoints,field,timestamp(),'x','z')  
        
    def plot2dPts(self,arr,npoints,field,tag,ax0,ax1):
#         fig = figure()
        i0 = ('x','y','z').index(ax0)
        i1 = ('x','y','z').index(ax1)
        scatter(arr[:npoints][field][:,i0],arr[:npoints][field][:,i1])
        xlabel('{}'.format(ax0))
        ylabel('{}'.format(ax1))
        name = '{}_{}_{}-{}'.format(tag,field,ax0,ax1)
        title(name)
        savefig(self.path+'/'+ name+'.pdf');
        os.system('cp {} {}'.format(self.path+'/'+ name+'.pdf',self.path+ '/latest{}-{}'.format(ax0,ax1)+'.pdf' ))
#         show()

