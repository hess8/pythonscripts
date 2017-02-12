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
11. Copy all the ind. points by trans and pt symmetry again. r
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
from scipy.spatial import Delaunay as delaunay, Voronoi as sci_voronoi
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
    reverseStructured,isInVoronoi)

def timestamp():
    return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())


class meshConstruct(): 
    '''Makes superperiodic structures and cuts'''
    from comMethods import readfile,writefile,trimSmall,areEqual,directFromCart,cartFromDirect
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
    from conjGradMin2 import (fmin_cg,minimize_cg,line_search_wolfe1,scalar_search_wolfe1)

    def __init__(self):
        self.strLen = 5
        self.indIn = []
        self.indOut = []
                 
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
        self.edgeFactor = 3.0
        self.rmin = 0.8*self.ravg #
        self.rmax = self.edgeFactor*self.ravg #cutoff for forces, neighbors in other cells. 
        print 'Typical spacing:', self.ravg
        self.initPoints() 
        self.initTets()
        self.movetoVoronoi()
        self.addBorders()
        self.delauPoints()
        self.fillTets()
        self.relax()
        self.plotPoints(self.points,self.npoints,'pos')
        sys.exit('stop')     
        return meshvecs, Nmesh, lattype, pfB, pf, status

    def movetoVoronoi(self):
        for i in range(self.npoints):
            self.points[i]['pos'] = intoVoronoi(self.points[i]['pos'],self.B)
    
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
        newPoint = zeros(1,dtype = [('dep', 'S8'),('inCell', 'bool'),
            ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('force', '3float')])
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
            self.indIn.append(ipoint)
        else:
            newPoint[0]['inCell'] = False
            inCell = False
            self.indOut.append(ipoint)
        self.points = concatenate((self.points,newPoint),axis = 0)
        self.sympoints = [str(pos[0])[:self.strLen]+str(pos[1])[:self.strLen]+str(pos[2])[:self.strLen]]
        return pos,inCell 

    def addPntsSym(self,pos,isponsor,ipoint,type = None):
        '''Add points by point-symmetry.  Check that each is far enough away 
        from previous ones.  If too close, choose the point midway, make this the
        the independent point, and run through all the sym operators again.
        
        If the sponsor point is in the cell, we will move symmetry partners into cell.
        If not, we will not map them.  
        '''
        newPoints = zeros(self.nops -1,dtype = [('dep', 'S8'),('inCell', 'bool'),
                ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('force', '3float')])
        doneSym = False
        nNew = 0
        while not doneSym: #symmetry dependents
            posList = [pos]
            ipos = 0
            for op in range(self.nops): 
                newpos = dot(self.symops[:,:,op],pos) #using voronoi cell. Sym ops don't take them out of cell. 
#                 if type == 'init':
#                     newpos = intoCell(dot(self.symops[:,:,op],pos),self.B)
#                 else: 
#                     newpos = dot(self.symops[:,:,op],pos)
                posList.append(newpos)
                [tooClose,closePos] = self.tooCloseList(posList,self.ravg/3.0)
                if not tooClose:
                    posStr = str(newpos[0])[:self.strLen]+str(newpos[1])[:self.strLen]+str(newpos[2])[:self.strLen]
                    if posStr not in self.sympoints:
                        self.sympoints.append(posStr)
                        newPoints[ipos]['pos'] = newpos
                        dposnew = self.directFromCart(self.B,newpos)    
                        newPoints[ipos]['dpos'] = dposnew
                        newPoints[ipos]['dep'] = op
                        newPoints[ipos]['sponsor'] = isponsor
                        df = dposnew.flatten()
                        inCell = (min(df) >= 0 and max(df) < 1 )
                        newPoints[ipos]['inCell'] = inCell
                        ipos += 1                          
                else:
                    #combine points, write over the last independent point, and start over on ops
                    print 'combine',newpos , closePos
                    pos = (newpos + closePos)/2.0
#                     if ipoint -1 in self.indIn: 
#                         self.indIn.remove(ipoint-1)
#                     else: 
#                         self.indOut.remove(ipoint-1)
                    self.points[ipoint-1]['pos'] = pos
                    dpos = self.directFromCart(self.B,pos)
                    self.points[ipoint-1]['dpos'] = dpos
                    self.points[ipoint-1]['inCell'] = inCell
                    self.sympoints = [str(pos[0])[:self.strLen]+str(pos[1])[:self.strLen]+str(pos[2])[:self.strLen]]
                    break                    
                    #keep dep the same: tet it came from...but this means that point may not be inside tet.
            doneSym = True
            for i in range(ipos):
                if newPoints[i]['inCell']: self.nCellPoints += 1                     
        self.points = concatenate((self.points,newPoints[:ipos]),axis = 0)
        return ipos
            


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
        '''Add points to the largest tets, until the 
        the desired number of points is obtained.  We add a point to the largest
        tetrahedron, then replicate it symmetrically, then resort tets by volume 
        (the Delaunay routine doesn't seem to make a triangulation that preserves
        symmetry.'''
        ipoint = self.npoints-1
        inewpt = 0
        if self.nops > 1:
            while (self.nCellPoints < self.nTarget): #add nops points at a time1:
                tet = self.tets[0]['tet'] #largest volume tet
                [pos0,inCell] = self.addPntTet(0,ipoint)
                ipoint += 1
                isponsor = ipoint
                nNew = self.addPntsSym(pos0,isponsor,ipoint,inCell)
                ipoint += nNew       
                self.npoints = ipoint + 1 
                print'plotting {} points'.format(self.npoints)
                self.plotPoints(self.points,'pos',str(ipoint+1),nNew,self.tets[0:1])
                #resort according to volume (new tets are now old), and then check if we need more points: 
                self.delauPoints() #find all tets and sort    
        else: # no symmetries
            while (self.nCellPoints < self.nTarget):
                itet = 0
                firstVol = self.tet[0]['vol'] 
                while self.tet[itet]['vol'] > firstVol/4.0:
                    self.addPntTet(itet, ipoint)
                    ipoint += 1
                self.delauPoints() #find all tets and sort   
                print'plotting {} points'.format(self.npoints)
                self.plotPoints(self.points,'pos',str(ipoint+1))   

    def initTets(self):
        '''Add points to the largest tets then replicate them symmetrically, 
        then resort tets by volume
        If any points are too close upon the symmetry operations, we combine them.  
        '''
        self.delauPoints()
        ipoint = self.npoints
        self.nCellPoints = 0 #don't count the VC vertices as points
        if self.nops > 1:
            while (ipoint < 8+8): #add nops points at a time1:
                tet = self.tets[0]['tet'] #largest volume tet
                [pos0,inCell] = self.addPntTet(0,ipoint)
                ipoint += 1
                isponsor = ipoint
                nNew = self.addPntsSym(pos0,isponsor,ipoint,'init') #We are working in original parallelpiped
                ipoint += nNew   
#                 print'plotting {} points'.format(ipoint)
#                 self.plotPoints(self.points,'pos',str(ipoint),self.tets[0:1])
                #resort according to volume (new tets are now old), and then check if we need more points: 
                self.delauPoints() #find all tets and sort    
        else: # no symmetries
            while (ipoint < 8+8):
                itet = 0
                firstVol = self.tets[0]['vol'] 
                while self.tets[itet]['vol'] > firstVol/4.0:
                    self.addPntTet(itet, ipoint)
                    ipoint += 1
                self.delauPoints() #find all tets and sort    
        
#         self.cellPoints = self.points[8:]  #  remove the original 8 tets, which were on the boundary of the cell.               
#         self.ncellPoints = len(self.cellPoints)
#         self.indIn = [i-8 for i in self.indIn]
        
#         self.cellPoints = self.points          
        self.nCellPoints = len(self.points)
        self.npoints = self.nCellPoints
        self.nIndIn = len(self.indIn)
        print'plotting {} points'.format(self.nCellPoints)
        self.plotPoints(self.points,'pos',str(self.nCellPoints),0,self.tets[:])
#         self.nCellPoints = ipoint + 1 - 8
        return       

    def initPoints(self):
        #self.subRandSym()
#         #use B cell vertices as the first 8 points
#         self.points = zeros(8,dtype = [('dep', 'S8'),('inCell', 'bool'),
#             ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('force', '3float')])
#         delta = 1e-6
#         a = 1-delta
#         z = 0+delta
#         for i,mult in enumerate([[z,z,z], [z,z,a],[z,a,z],[a,z,z],[z,a,a],[a,z,a],[a,a,z],[a,a,a]]):
#             dpos = array(mult)
#             self.points[i]['dpos'] = dpos
#             self.points[i]['pos'] = dot(self.B,dpos)
#             self.points[i]['inCell'] = True
#             #will remove these 8 points later, so don't worry about the rest of the info
        
        #make a voronoi cell at the origin.  Use its vertices as the initial points
        Lvecs = []
        for i in [0,-1,1]:
            for j in [0,-1,1]:
                for k in [0,-1,1]:
                    Lvecs.append(dot(self.B,array([i,j,k])))
        vorPoints = sci_voronoi(array(Lvecs)).vertices
        nvp = len(vorPoints)
        print 'number of vornoi vertices', nvp    
        self.points = zeros(nvp+1,dtype = [('dep', 'S8'),('inCell', 'bool'),
            ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('force', '3float')])
        
        self.points[0]['pos'] = array([0.0,0.0,0.0])
        ipos = 1
        for i in range(nvp):
            self.points[i]['pos'] = vorPoints[i]
            self.points[i]['inCell'] = True
            ipos += 1
        self.npoints = len(self.points)
        print'plotting {} points'.format(self.npoints)
        self.plotPoints(self.points,'pos',str(self.npoints))

        return

    def addBorders(self):
        '''Add bordering points by translation symmetry to the Voronoi cell,
        keeping only those within an edge distance from the borders.
        We make a translational copy of the VC by the 27-1 translational vectors R
        nearest the center.  u_R is a unit vector in the direction of R.
        A new points has position r.  if dot(r,u_R)> R/2, 
        then it's outside the VC.  So we can use the conditon  dot(r,u_R) <  R/2 + edge
        to choose those that are in the border region.            
         '''
        ipoint = self.npoints
        for i in [0,-1,1]:
            for j in [0,-1,1]:
                for k in [0,-1,1]:
                    for ipos in range(self.nCellPoints):
                        if not (i==0 and j==0 and k ==0):
                           transVec = dot(self.B,array([i,j,k]))
                           normTV = norm(transVec)
                           unitVec = transVec/normTV
                           newPos = self.points[ipos]['pos'] + transVec
                           if dot(newPos,unitVec) <= normTV/2.0 + self.rmax: #rmax gives the border thickness   
                               newPoint = zeros(1,dtype = [('dep', 'S8'),('inCell', 'bool'),
                                    ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('force', '3float')])
                               newPoint['dep'] = '{},{},{}'.format(i,j,k) # this point is dependent by translation
                               newPoint['sponsor'] = ipos #The sponsor may itself be dependent.  
                               newPoint['inCell'] = False                        
                               newPoint['pos'] = newPos
#                                self.points[ipoint]['dpos'] = newDirPos
                               self.points = concatenate((self.points,newPoint),axis = 0)
                               ipoint += 1
        self.npoints = ipoint
        self.points = delete(self.points,s_[self.npoints:],0)
        if len(self.points) != self.npoints: sys.exit('Stop.  Error in numbering points in expanded cell')
        print 'Init points plus border: {}'.format(self.npoints)
        self.plotPoints(self.points,'pos','expanded')
        return   
                    

        
#     def combineNear(self): 
#         '''Doesn't preserve symmetry, so not using it yet'''
#         self.cellPoints2 = zeros(self.nCellPoints+self.nops,dtype = [('dep', 'S8'),('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),('cell', '3int')])
# #         combine = []
#         for jpos in range(self.nCellPoints):
#             for kpos in range(jpos+1,self.nCellPoints):
#                 r = norm(self.cellPoints[jpos]['pos'] - self.cellPoints[kpos]['pos']) < self.rmin
#                 if 0.0 < norm(self.cellPoints[jpos]['pos'] - self.cellPoints[kpos]['pos']) < self.rmin:
#                     print jpos,kpos
#                     print 'r', norm(self.cellPoints[jpos]['pos'] - self.cellPoints[kpos]['pos'])
#                     if kpos in self.labels: self.labels.remove(kpos)
#                     if kpos in self.indIn: self.indIn.remove(kpos)
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
        self.cellPoints = zeros(self.nTarget+self.nops,dtype = [('dep', 'S8'),
            ('sponsor', 'int8'),('pos', '3float'),('dpos', '3float'),
            ('force', '3float')])
        self.directPos = zeros((self.nTarget,3))
        self.subrand = rand(3)
        afact = sqrt(2)
        bfact = sqrt(3)
        cfact = sqrt(5)
        self.indIn = []
        self.indOut = []
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
            if ipos > 0 and self.tooCloseInit(pos,dpos,ipos,edges):
                print 'Independent point {} was too close. Trying again'.format(ipos)
                ipos = itemp
                continue #back to while
            self.cellPoints[ipos]['dep'] = 'I' #independent
            self.cellPoints[ipos]['sponsor'] = ipos                              
            self.cellPoints[ipos]['pos'] = pos 
            self.cellPoints[ipos]['dpos'] = dpos #dpos will always be the direct coordinates in the ORIGINAL cell
#             self.labels.append(ipos)
            self.indIn.append(ipos)
            iInd = ipos
            #now all symmetry operations will keep the point in the BZ 
            sympoints = [str(pos[0])[:self.strLen]+str(pos[1])[:self.strLen]+str(pos[2])[:self.strLen]]
            for op in range(self.nops):
                newpos = dot(self.symops[:,:,op],pos)
                posStr = str(newpos[0])[:self.strLen]+str(newpos[1])[:self.strLen]+str(newpos[2])[:self.strLen]
                if posStr not in sympoints:
                    dposnew = self.directFromCart(self.B,newpos) 
                    if self.tooCloseInit(newpos,dposnew,ipos,edges):
                        print 'Dependent point {} was too close. Trying again'.format(ipos)
                        ipos = itemp
                        self.indIn.remove(iInd)
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
        self.cellPoints = delete(self.cellPoints,s_[self.nCellPoints:],0) #removing unused rows
        if len(self.cellPoints) != self.nCellPoints: sys.exit('Stop.  Error in numbering cellPoints')
        self.nIndIn = len(self.indIn) 
        self.plotPoints(self.cellPoints,'pos','voronoi')
#         self.plotPoints(self.cellPoints,'dpos')
        print 'Points in unit cell:',self.nCellPoints   
        print 'Independent points:',self.nIndIn                         
        return

    def tooCloseList(self,posList,tol):
        '''Checks the last pos vs the others'''
        for partner in posList[:len(posList)-1]:
            r = norm(posList[-1]-partner)
            if 1e-6 < r < tol:
                print 'Too close to another:',r
                return True, partner
        return False,[] 
    
    def tooCloseInit(self,pos,dpos,ipos,edges):
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
        for i in self.indIn:
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
                plot([pair[0][i0],pair[1][i0]],[pair[0][i1],pair[1][i1]]) #draw line for tet edge
        xlabel('{}'.format(ax0))
        ylabel('{}'.format(ax1))
        name = '{}_{}_{}-{}'.format(tag,field,ax0,ax1)
        title(name)
        fig.savefig(self.path+'/'+ name+'.pdf');
        os.system('cp {} {}'.format(self.path+'/'+ name+'.pdf',self.path+ '/latest{}-{}'.format(ax0,ax1)+'.pdf' ))
        close()

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
        self.oldindVecs = [self.points[i]['pos'] for i in self.indIn]
        self.indInComps = array([self.points[i]['pos'] for i in self.indIn]).flatten()
#         initialGuess = [i+ 0.5*self.rmin*rand() for i in self.indInComps]  #just add noise to initial positions
        initialGuess = self.indInComps
        self.count = 0
        
        epsilon = self.ravg/100
        out = self.fmin_cg(initialGuess,epsilon)
        print out
#         for i in range(27):
#             print i ,self.indInComps[i]
#         vectors = self.indInComps.reshape((self.npoints,3))
#         for i in range(9):
#             print i ,vectors[i,:]        
        return
    
    def updatePoints(self,indVecs):
        for ipoint in range(self.npoints):
            if ipoint in self.indIn:
                place = self.indIn.index(ipoint)
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
                    sys.exit('Stop. Point {} is labeled as independent, but is not in self.indIn'.format(ipoint))
                elif ',' not in dep: # symmetry dependent
                    self.points[ipoint]['pos'] = intoCell(dot(self.symops[:,:,int(dep)],self.points[sponsor]['pos']),self.B)
                else: #translation
                    ms = array([int(i) for i in dep.split(',')])
                    svec = self.points[sponsor]['pos']
                    self.points[ipoint]['pos'] = svec + dot(self.B,ms)       
 
    def enerGrad(self,indComps):
        '''restrict the energy sum to the pairs that contain independent points'''
#         print 'oldindvecs',self.oldindVecs
        indVecs = indComps.reshape((self.nIndIn,3))
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
        for ipos in range(self.nIndIn):
            force = zeros(3)
            rs = zeros(self.npoints,dtype = [('r', 'float'),('rij', '3float'),('force', '3float')])
                      
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
            
            if ipos in self.indIn:
                rs =sort(rs, order='r') 
#                 print ipos, rs[:20]['r'] 
#                 print rs[:20]['rij']
                print ipos, 'pos', self.points[ipos]['pos']
                print 'nerst',self.points[nearest]['pos'],nearest,rnearest 
                print 'vectr', trimSmall(self.points[nearest]['pos'] - self.points[ipos]['pos'])
                print 'force', ipos, self.points[ipos]['force']    
        ener = ener/self.nIndIn
        grad = -self.points[:self.nIndIn]['force'].flatten()/self.nIndIn 
        print 'energy:',self.count, ener
        print
#         if self.count == 19:
#             self.plotPoints(self.points,'pos')
        self.count += 1
        #now update all the dependent positions
      #  ! Note:  need to update the positions at each step!  Do we have to get inside
#         return ener#         
        return ener, grad
