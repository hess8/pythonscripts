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
from numpy import (mod,dot,cross,transpose, rint,floor,ceil,zeros,array,sqrt,
                   average,std,amax,amin,int32,sort,count_nonzero,
                   delete,mean,square,argmax,argmin,insert,s_,concatenate)
# from scipy.optimize import fmin_cg
from scipy.spatial import Delaunay as delaunay
from numpy.linalg import inv, norm, det
from numpy.random import rand
from copy import copy,deepcopy
from sched import scheduler
from itertools import chain
from matplotlib.pyplot import (subplots,savefig,imshow,close,plot,title,xlabel,
                               ylabel,figure,show,scatter,triplot)
import matplotlib.image as mpimg
import datetime
from _ast import operator
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
sys.path.append('/fslhome/bch/graphener/graphener')


# from conjGradMin.optimize import fmin_cg

from kmeshroutines import (svmesh,svmeshNoCheck,svmesh1freedir, lattice_vecs, lattice, surfvol,
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double,
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection,intoVoronoi,intoCell,
    reverseStructured)

def timestamp():
    return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())


class meshConstruct(): 
    '''Makes superperiodic structures and cuts'''
    from comMethods import readfile,writefile,trimSmall,areEqual,directFromCart,cartFromDirect
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
    from conjGradMin2 import (fmin_cg,minimize_cg,line_search_wolfe1,scalar_search_wolfe1)

    def __init__(self):
        '''init'''
         
    def relaxMeshSym(self,B,targetNmesh,path):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        self.B = B
        print 'B',B
        [self.symops,self.nops] = getGroup(self.B)
        self.nTarget = targetNmesh
        self.path = path
        print 'Number of desired points:', targetNmesh
        print 'Symmetry operations:', self.nops
        vol = abs(det(B))
        self.ravg = (vol/targetNmesh)**(1/3.0)
        self.edgeFactor = 1.5
        self.rmin = 0.8*self.ravg #
        self.rmax = self.edgeFactor*self.ravg #cutoff for forces, neighbors in other cells. 
        print 'Typical spacing:', self.ravg
        self.initPoints() 
        self.transPoints()
        self.movetoVoronoi()
        self.delauPoints()
        self.fillTets()
        self.relax()
        self.plotPoints(self.points,self.npoints,'pos')
        sys.exit('stop')     
        return meshvecs, Nmesh, lattype, pfB, pf, status

    def movetoVoronoi(self):
        for i in self.npoints:
            self.points['pos'] = self.toVoronoi(self.points['pos'])
    
    def volTet(self,vertPoints): 
        return dot((vertPoints[0]-vertPoints[3]),cross((vertPoints[1]-vertPoints[3]),(vertPoints[2]-vertPoints[3])))/6.0
    
    def delauPoints(self):
#         print 'len self.points', len(self.points[:self.npoints]['pos']) 
        tri = delaunay(self.points['pos'])
        self.ntets = tri.nsimplex
        self.tets = zeros(self.ntets, dtype = [('tet', '4int'),('vol', 'float'),('anyInCell', 'bool')])
        self.tets['tet'] = tri.simplices #vertex labels
        self.tets['anyInCell'] = False
        for itet, tet in enumerate(self.tets[:]['tet']):
            for ivert in range(4): #find those with at least one vertex in the cell
                if self.points[tet[ivert]]['inCell']:
                    self.tets[itet]['anyInCell'] = True
                    break
            if self.tets[itet]['anyInCell']:
                self.tets[itet]['vol'] = abs(self.volTet([self.points[iv]['pos'] for iv in tet]))
        self.tets.sort(order='vol');self.tets = reverseStructured(self.tets)

    def addPntTet(self,itet,ipoint):
        tet = self.tets[itet]['tet'] 
        vp = array([self.points[iv]['pos'] for iv in tet])
        vd = array([self.points[iv]['dpos'] for iv in tet])
        newPoint = zeros(1,dtype = [('label', 'int8'),('dep', 'S8'),('inCell', 'bool'),
            ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('force', '3float')])
        newPoint[0]['label'] = ipoint
        pos = sum(vp)/4.0 #center of mass
        newPoint[0]['pos'] = pos
#         print 'newpos',newPoint[0]['pos']
        dpos = sum(vd)/4.0 
        newPoint[0]['dpos'] = dpos
        newPoint[0]['dep'] = 'tet'.format(tet[0],tet[1],tet[2],tet[3]) 
        if min(dpos.flatten()) >= 0 and max(dpos.flatten()) < 1.0:
            newPoint[0]['inCell'] = True
            inCell = True
            self.nCellPoints += 1
        else:
            newPoint[0]['inCell'] = False
            inCell = False
        self.points = concatenate((self.points,newPoint),axis = 0)
        self.sympoints = [str(pos[0])[:self.strLen]+str(pos[1])[:self.strLen]+str(pos[2])[:self.strLen]]
        return pos,inCell 

    def addPntSym(self,pos0,op,isponsor,ipoint):
        '''In Voronoi cell'''
        newpos = dot(self.symops[:,:,op],pos)  
        posStr = str(newpos[0])[:self.strLen]+str(newpos[1])[:self.strLen]+str(newpos[2])[:self.strLen]
        sympoints = [str(pos[0])[:self.strLen]+str(pos[1])[:self.strLen]+str(pos[2])[:self.strLen]]
        if posStr not in self.sympoints:
            newPoint = zeros(1,dtype = [('label', 'int8'),('dep', 'S8'),('inCell', 'bool'),
                ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('force', '3float')])
            sypoints.append(sympoints)
            newPoint[0]['pos'] = newpos
            dposnew = self.directFromCart(newpos)
            newPoint[0]['dpos'] = dposnew
            newPoint[0]['dep'] = op
            newPoint[0]['sponsor'] = isponsor
            newPoint[0]['inCell'] = True
            if min(dposnew.flatten()) >= 0 and max(dposnew.flatten()) < 1.0:
                newPoint[0]['inCell'] = True
                self.nCellPoints += 1
                self.ind.append(ipoint)
            else:
                newPoint[0]['inCell'] = False
            self.points = concatenate((self.points,newPoint),axis = 0)
            return True
        else:
            return False

#     def addTets(self,itet):
#         newTets = zeros(4,dtype = [('tet', '4int'),('vol', 'float'),('anyInCell', 'bool')])
#         for i in range(4):
#             newTets[0]['tet'] = [ipoint,tet[0],tet[1],tet[2]]
#             newTets[1]['tet'] = [ipoint,tet[0],tet[1],tet[3]]
#             newTets[2]['tet'] = [ipoint,tet[1],tet[2],tet[3]]
#             newTets[3]['tet'] = [ipoint,tet[2],tet[3],tet[1]]
#             newTets[0]['anyInCell'] = (True in [inCell,self.points[0]['inCell'],self.points[1]['inCell'],self.points[2]['inCell']])
#             newTets[0]['anyInCell'] = (True in [inCell,self.points[0]['inCell'],self.points[1]['inCell'],self.points[3]['inCell']])
#             newTets[0]['anyInCell'] = (True in [inCell,self.points[1]['inCell'],self.points[2]['inCell'],self.points[3]['inCell']])
#             newTets[0]['anyInCell'] = (True in [inCell,self.points[2]['inCell'],self.points[3]['inCell'],self.points[1]['inCell']])
#         self.tets = concatenate((self.tets,newTets),axis = 0)
#         self.ntets += 4 #total number (next index of tet to add)
    
    def fillTets(self):
        '''Add points to the largest tet, which creates 4 new ones, until the 
        the desired number of points is obtained.  We add a point to the largest
        tetrahedron, then replicate it symmetrically, then resort tets by volume 
        (the Delaunay routine doesn't seem to make a triangulation that preserves
        symmetry.'''
        ipoint = self.npoints-1
        ivol = 0 #may be used to space plotting
        inewpt = 0
        if self.nops > 1:
            while (self.nCellPoints < self.nTarget): #add nops points at a time1:
                tet = self.tets['tets'][0] #largest volume tet
                pos0 = self.addPntTet(itet,ipoint)
                ipoint += 1
                for op in range(self.nops): #symmetry dependents
                    if self.addPntSym(pos0,op,ipoint):
                        ipoint += 1
                self.npoints = self.ipoint + 1 
                print'plotting {} points'.format(self.npoints)
                self.plotPoints(self.points,'pos',str(ipoint+1),self.tets[0:1])
                #resort according to volume (new tets are now old), and then check if we need more points: 
                self.delauPoints() #find all tets and sort    
        else: # no symmetries
            while (self.nCellPoints < self.nTarget):
                itet = 0
                firstVol = self.tet[itet]['vol'] 
                while self.tet[itet]['vol'] > firstVol/4.0:
                    self.addPntTet(itet, ipoint)
                    ipoint += 1
                self.delauPoints() #find all tets and sort   
                print'plotting {} points'.format(self.npoints)
                self.plotPoints(self.points,'pos',str(ipoint+1))                    

    def relax(self):
        '''Conjugate gradient relaxation of the potential energy.
        See www.idi.ntnu.no/~elster/tdt24/tdt24-f09/cg.pdf 
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
        self.oldindVecs = [self.points[i]['pos'] for i in self.ind]
        self.indComps = array([self.points[i]['pos'] for i in self.ind]).flatten()
#         initialGuess = [i+ 0.5*self.rmin*rand() for i in self.indComps]  #just add noise to initial positions
        initialGuess = self.indComps
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
                self.points[ipoint]['pos'] = intoCell(indVecs[place,:],self.B) 
#                 self.points[ipoint]['pos'] = indVecs[place,:]          
            else:
                dep = self.points[ipoint]['dep']
                sponsor = self.points[ipoint]['sponsor']
#                 print 'test',ipoint
#                 print 'dep' , dep
#                 print 'spons',sponsor
#                 print
                if dep == 'I':
                    sys.exit('Stop. Point {} is labeled as independent, but is not in self.ind'.format(ipoint))
                elif ',' not in dep: # symmetry dependent
                    self.points[ipoint]['pos'] = intoCell(dot(self.symops[:,:,int(dep)],self.points[sponsor]['pos']),self.B)
                else: #translation
                    ms = array([int(i) for i in dep.split(',')])
                    svec = self.points[sponsor]['pos']
                    self.points[ipoint]['pos'] = svec + dot(self.B,ms)       
 
    def enerGrad(self,indComps):
        '''restrict the energy sum to the pairs that contain independent points'''
#         print 'oldindvecs',self.oldindVecs
        indVecs = indComps.reshape((self.nInd,3))
        self.oldindVecs = indVecs
        #update all positions by symmetry and translation
        self.oldPoints = deepcopy(self.points)
        self.updatePoints(indVecs)
        #indVecs now cannot be used for current points because they might be out of cell.
#         print 'movement'
#         for i in range(self.npoints):
# #             print 'old',i,self.oldPoints[i]['pos']
# #             print 'new',i,self.points[i]['pos']
# #             print
#             move = norm(self.points[i]['pos']-self.oldPoints[i]['pos'])
#             if move > 1e-6 :
#                 print i,move
        
        ener = 0.0
        self.power = 2.0
#         scale = 1
        for ipos in range(self.nInd):
            force = zeros(3)
            rs = zeros(self.npoints,dtype = [('label', 'int8'),('r', 'float'),('rij', '3float'),('force', '3float')])
                      
            nearest = -1
            rnearest = 100.0
            for jpos in range(self.npoints):
                rij = self.points[jpos]['pos'] - self.points[ipos]['pos']
                r = norm([rij])
                if r > 1e-4*self.rmin:
                    if r<rnearest:
                        nearest = jpos
                        rnearest = r
#                    
#                     epair = r**2 #attractive springs...dumb because those farthest away influence most
#                     force += rij #attractive springs)
                    
                    epair = (self.ravg/r)**self.power
                    forcepair = -self.power*epair/r * rij/r
                    force += forcepair
#                     print 'jpos',jpos,r,forcepair
                    ener += epair
                    rs[jpos]['r'] = r
                    rs[jpos]['rij'] = rij
                    rs[jpos]['force'] = forcepair
              
                    
            self.points[ipos]['force'] = trimSmall(force)
            
            if ipos in self.ind:
                rs =sort(rs, order='r') 
#                 print ipos, rs[:20]['r'] 
#                 print rs[:20]['rij']
                print ipos, 'pos', self.points[ipos]['pos']
                print 'nerst',self.points[nearest]['pos'],nearest,rnearest 
                print 'vectr', trimSmall(self.points[nearest]['pos'] - self.points[ipos]['pos'])
                print 'force', ipos, self.points[ipos]['force']    
        ener = ener/self.nInd
        grad = -self.points[:self.nInd]['force'].flatten()/self.nInd 
        print 'energy:',self.count, ener
        print
#         if self.count == 19:
#             self.plotPoints(self.points,'pos')
        self.count += 1
        #now update all the dependent positions
      #  ! Note:  need to update the positions at each step!  Do we have to get inside
#         return ener#         
        return ener, grad

    def transPoints(self):
        '''These points are in the parallelpiped ORIGINAL cell.
        Make copies of points in 3x3x3 supercell, but keep only those that 
        are within rc of edges'''
        self.points = zeros(self.nCellPoints*27,dtype = [('label', 'int8'),('dep', 'S8'),('inCell', 'bool'),('sponsor', 'int8'),('pos', '3float'),
                                                     ('dpos', '3float'),('force', '3float')])
        aEdge = self.rmax/norm(self.B[:,0]); bEdge = self.rmax/norm(self.B[:,1]); cEdge = self.rmax/norm(self.B[:,2])
        ipoint = 0
        for i in [0,-1,1]:
            for j in [0,-1,1]:
                for k in [0,-1,1]:
                    for ipos in range(self.nCellPoints):
                        dpos = self.cellPoints[ipos]['dpos'] #dpos will always be the direct coordinates in the original cell
                        #move to parallepiped cell
                        pos = dot(self.B,dpos)
                        self.cellPoints[ipos]['pos']  = pos #was voronoi position before                
                        newDirPos = dpos + array([i,j,k])
                        if  (-aEdge < newDirPos[0] < 1 + aEdge) and (-bEdge < newDirPos[1] < 1 + bEdge) and (-cEdge < newDirPos[2] < 1 + cEdge):
                            newPos = dot(self.B,newDirPos)
#                             print 'ipoint,ipos',ipoint,ipos, newPos
#                             print '\tdpos',dpos
#                             self.points[ipoint]['label'] = ipoint
#                             self.labels.append(ipoint)
                            if i == 0 and j == 0 and k == 0: #original cell
                                self.points[ipoint]['dep'] = self.cellPoints[ipoint]['dep'] 
                                self.points[ipoint]['sponsor'] = self.cellPoints[ipoint]['sponsor']
                                self.points[ipoint]['inCell'] = True
                            else: #then this point is dependent by translation
                                self.points[ipoint]['dep'] = '{},{},{}'.format(i,j,k)
                                self.points[ipoint]['sponsor'] = ipos #The sponsor may itself be dependent.  
                                self.points[ipoint]['inCell'] = False                        
                            self.points[ipoint]['pos'] = newPos
                            self.points[ipoint]['dpos'] = newDirPos
                            ipoint += 1
        self.npoints = ipoint
        self.points = delete(self.points,s_[self.npoints:],0)
        if len(self.points) != self.npoints: sys.exit('Stop.  Error in numbering points in expanded cell')
        print 'All points: {}'.format(self.npoints)
        self.plotPoints(self.points,'pos','expanded')
#         self.plotPoints(self.points,'dpos')
        return   
    
    def initPoints(self):
        self.subRandSym()  

        
#     def combineNear(self): 
#         '''Doesn't preserve symmetry, so not using it yet'''
#         self.cellPoints2 = zeros(self.nCellPoints+self.nops,dtype = [('label', 'int8'),('dep', 'S8'),('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('cell', '3int')])
# #         combine = []
#         for jpos in range(self.nCellPoints):
#             for kpos in range(jpos+1,self.nCellPoints):
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
#         self.nCellPoints = ipos
#         self.cellPoints[label] = self.cellPoints2[ipos]
#         print 'Points generated:',self.nCellPoints
#         return
                      
    def subRandSym(self):
        '''These points are in the VORONOI cell'''
#         self.labels = [] # list of active points labels
        self.cellPoints = zeros(self.nTarget+self.nops,dtype = [('label', 'int8'),('dep', 'S8'),
            ('sponsor', 'int8'),('pos', '3float'),('vpos', '3float'),('dpos', '3float'),
            ('force', '3float')])
        self.directPos = zeros((self.nTarget,3))
        self.subrand = rand(3)
        afact = sqrt(2)
        bfact = sqrt(3)
        cfact = sqrt(5)
        self.ind = []
        ipos = -1
        edgeFactor = 10
        aEdge = self.rmax/norm(self.B[:,0])/edgeFactor; bEdge = self.rmax/norm(self.B[:,1])/edgeFactor; cEdge = self.rmax/norm(self.B[:,2])/edgeFactor        
        edges = [aEdge,bEdge,cEdge]
#         while ipos < self.nTarget: #do it once for now.
        while ipos < 0:
            itemp = ipos
            ipos += 1
            self.subrand = mod(self.subrand + array([afact,bfact,cfact]), array([1,1,1]))
            dpos = self.subrand
            pos = dot(self.B,dpos)
            pos = intoVoronoi(pos,self.B)#now in Voronoi cell
            print 'ind point', ipos,pos
            if ipos > 0 and self.tooClose(pos,dpos,ipos,edges):
                print 'Independent point {} was too close. Trying again'.format(ipos)
                ipos = itemp
                continue #back to while
            self.cellPoints[ipos]['dep'] = 'I' #independent
            self.cellPoints[ipos]['sponsor'] = ipos                              
            self.cellPoints[ipos]['pos'] = pos 
            self.cellPoints[ipos]['dpos'] = dpos #dpos will always be the direct coordinates in the ORIGINAL cell
#             self.labels.append(ipos)
            self.ind.append(ipos)
            iInd = ipos
            #now all symmetry operations will keep the point in the BZ 
            self.strLen = 5
            sympoints = [str(pos[0])[:self.strLen]+str(pos[1])[:self.strLen]+str(pos[2])[:self.strLen]]
            for op in range(self.nops):
                newpos = dot(self.symops[:,:,op],pos)
                posStr = str(newpos[0])[:self.strLen]+str(newpos[1])[:self.strLen]+str(newpos[2])[:self.strLen]
                if posStr not in sympoints:
                    dposnew = self.directFromCart(self.B,intoCell(newpos,self.B)) 
                    if self.tooClose(newpos,dposnew,ipos,edges):
                        print 'Dependent point {} was too close. Trying again'.format(ipos)
                        ipos = itemp
                        self.ind.remove(iInd)
                        break
                    sympoints.append(posStr)
                    ipos += 1
#                     self.labels.append(ipos)
#                     self.dep.append(ipos)
                    self.cellPoints[ipos]['dep'] = '{}'.format(str(op))
                    self.cellPoints[ipos]['sponsor'] = iInd 
#                     newpos=intoVoronoi(newpos,self.B)                                
                    self.cellPoints[ipos]['dpos'] = dposnew
#                     checkpos = intoVoronoi(newpos,self.B)
#                     if norm(checkpos-newpos)> 1e-6:
#                         print 'New Voronoi cell pos for',newpos,checkpos                   
                    self.cellPoints[ipos]['pos'] = newpos  
#                     print 'ipos',ipos, newpos 
#                     print '\tdpos',self.cellPoints[ipos]['dpos'] 
#                 if ipos == itemp: #some point is too close...try another independent point
#                     break                
        self.nCellPoints = ipos + 1
        self.cellPoints = delete(self.cellPoints,s_[self.nCellPoints:],0)
        if len(self.cellPoints) != self.nCellPoints: sys.exit('Stop.  Error in numbering cellPoints')
        self.nInd = len(self.ind) 
        self.plotPoints(self.cellPoints,'pos','voronoi')
#         self.plotPoints(self.cellPoints,'dpos')
        print 'Points in unit cell:',self.nCellPoints   
        print 'Independent points:',self.nInd                         
        return
    
    def tooClose(self,pos,dpos,ipos,edges):
        for i in range(3):
            if abs(dpos[i]-0.5)> 0.5-edges[i]:
                print 'Too close to cell boundary'
                return True
        for testpos in self.cellPoints[:ipos]['pos']:
            if norm(pos-testpos) < self.ravg/5:
                print 'Too close to another:',norm(pos-testpos), testpos
                return True
        return False
    
    def plotPoints(self,pts,field = 'pos',tag = timestamp(),highlight = 0,tets = []):
        self.plot2dPts(pts,field,tag,'x','y',highlight,tets) 
        self.plot2dPts(pts,field,tag,'x','z',highlight,tets)  
        
    def plot2dPts(self,pts,field,tag,ax0,ax1,highlight,tets):
        fig = figure()
        i0 = ('x','y','z').index(ax0)
        i1 = ('x','y','z').index(ax1)
        scatter(pts[:][field][:,i0],pts[:][field][:,i1])
        if highlight > 0:
            scatter(pts[-highlight:][field][:,i0],pts[-highlight:][field][:,i1],color = 'g')
        ind0 = []
        ind1 = []
        for i in self.ind:
            ind0.append(pts[i][field][i0])
            ind1.append(pts[i][field][i1])            
        scatter(ind0,ind1,color = 'r') 
        if len(tets)>0: 
            pairs = []
            for itet, tet in enumerate(tets[:]['tet']):
                pairs.append([pts[tet[0]][field],pts[tet[1]][field]])
                pairs.append([pts[tet[0]][field],pts[tet[2]][field]])
                pairs.append([pts[tet[0]][field],pts[tet[3]][field]])
                pairs.append([pts[tet[1]][field],pts[tet[2]][field]])
                pairs.append([pts[tet[1]][field],pts[tet[3]][field]])
                pairs.append([pts[tet[2]][field],pts[tet[3]][field]])
            for pair in pairs:
                plot([pair[0][i0],pair[1][i0]],[pair[0][i1],pair[1][i1]],'k-') #draw line for tet edge
        xlabel('{}'.format(ax0))
        ylabel('{}'.format(ax1))
        name = '{}_{}_{}-{}'.format(tag,field,ax0,ax1)
        title(name)
        fig.savefig(self.path+'/'+ name+'.pdf');
        os.system('cp {} {}'.format(self.path+'/'+ name+'.pdf',self.path+ '/latest{}-{}'.format(ax0,ax1)+'.pdf' ))
        close()

