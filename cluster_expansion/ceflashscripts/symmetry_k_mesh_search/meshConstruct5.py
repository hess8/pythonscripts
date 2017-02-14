'''
1. Find the  NL point symmetries of the lattice 
2. Create the Latice iBZ: liBZ
3. Find Voronoi BZ
# 3.5. use any mirror planes to slice the VBZ
4. Populate the BZ with Nops copies. 
5. Find the NC symmetries of the crystal
6. keep NC copies that are contiguous: LBZ

All matrices store vectors as COLUMNS
'''


# '''
# Builds a mesh from the symmetry operations and random starting points.
# There is no lattice
# 0. The mesh should be optimized as an atomic space superlattice (transform of k-space)
# 0.1 Choose a 
# 1. nCoarse sub-random points in the BZ cell parallelpiped.  
#     We don't need to work in Voronoi cell version because we are running in 3x3x3 supercell of the BZ
#     https://en.wikipedia.org/wiki/Low-discrepancy_sequence#Additive_recurrence for subRandSym 
# 2. Define cutoff as rc = (nCoarse/Vcell)^1/3
# 3. Make all copies due to translation and symmetry operations
# 4. In neighboring cells (3x3x3 construction of plpd), keep only those that 
# are within e.g. 1.5 * rc of the cell boundaries
# 5. If two points are less than e.g. 0.8*rc away from each other, combine them 
# into a point at their center. 
# 5.5 if we have too many delete those that are closest to their neighbors
# 5.6 If we have too few, get the Delaunay tetrahedra, and add to the center 
# of the largest tets. 
# 6. Keep track of independent points (all are "sponsors") and symmetry partners
# Dependent points (which can also be sponsors) have "dep" label:  if by translations: "-1,0,1". 
# if by symmetry, a single vecitive integer "3" that gives the operation
# 7. Calculate the force on each ind. point from all the points that lie within
# a distance rForce.  Use a power law: Fi = k*Sum(r_vec_i,j/ rij^p
# 8.  Move the ind. points in a vector proportional to the force
# 9.  Copy all the ind. points by trans and pt symmetry again. 
# 10. repeat from 7 until the force on each point is less than minForce
# 11. Copy all the ind. points by trans and pt symmetry again. r
# ________
# Fine tetrahedral mesh
# ________
# 1. Find Delaunay tetrahedra
# 2. Subdivide each tet that has an independent point in it with a formula that 
# minimizes the S/V of the 6 new tets. Call these independent children, 
# even though they might have parents that are not independent
# 3. Copy these by trans and point sym, removing duplicates
# 4. Subdivide the indChildren (repeat 2 and 3) until we have the number of points we want.
# 5. By symmetry, find the number of truly independent among all points, and their weights. 
# 6. Write these to KPOINTS: https://cms.mpi.univie.ac.at/vasp/vasp/Entering_all_k_points_explicitly.html 
# until we have the number of points we want. 
# 
# All matrices store vectors as COULUMNS
# '''
import os, subprocess,sys,re,time
from numpy import (mod,dot,cross,transpose, rint,floor,ceil,zeros,array,sqrt,
                   average,std,amax,amin,int32,sort,count_nonzero,arctan2,
                   delete,mean,square,argmax,argmin,insert,s_,concatenate,all,
                   trace,where,real,allclose,sign,pi,imag,identity)
# from scipy.optimize import fmin_cg
from scipy.spatial import Delaunay as delaunay, Voronoi as sci_voronoi, ConvexHull as convexH
from numpy.linalg import inv, norm, det, eig
from numpy.random import rand
from copy import copy,deepcopy
from sched import scheduler
from itertools import chain, combinations, permutations
from matplotlib.pyplot import (subplots,savefig,imshow,close,plot,title,xlabel,
                               ylabel,figure,show,scatter,triplot)
import matplotlib.image as mpimg
import datetime
from _ast import operator
from pip._vendor.html5lib.constants import rcdataElements
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts')
sys.path.append('/fslhome/bch/graphener/graphener')

from meshpy.tet import MeshInfo, build, Options


# from conjGradMin.optimize import fmin_cg

from kmeshroutines import (svmesh,svmeshNoCheck,svmesh1freedir, lattice_vecs, lattice, surfvol,
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double,
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection,intoVoronoi,intoCell,
    reverseStructured,isInVoronoi,areParallel, addVec)

def timestamp():
    return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())

def threePlaneIntersect(rRows):  
    '''This routine is for three planes that will intersect at only one point, 
    as borders of facets, as we have for the Voronoi cell boundaries. 
    Planes are given by normals to ro, at the point ro.  All points r in the
    plane obey dot(r,ro) = ro^2 = dot(ro,ro)
      If they intersect, then 
        inv[[[xo,yo,zo]
        [x1,y1,z1]
        [x2,y2,z3]] }  (ro^2,r1^2,r2^2) has a solution
    '''
#     rRows = array(vecs3)
    rSq = array([dot(rRows[0],rRows[0]),dot(rRows[1],rRows[1]),dot(rRows[2],rRows[2])])
    try:
        invrRows = inv(rRows)
        invOK = True
    except:
        invOK = False
        return None
    if invOK:
        point = trimSmall(dot(invrRows,rSq))  
        if norm(point) < 100:
            return point
        else:
            return None
        
def orderAngle(facet,labels):
        '''get the angle that each vector is, in the plane of the facet.'''
        rcenter = sum(facet)/len(facet)
        xunitv =  (facet[0] - rcenter)/norm(facet[0] - rcenter)
#         yunitv = cross(xunitv, rcenter)/norm(cross(xunitv, rcenter))    
        vec = (facet[1] - rcenter) - dot((facet[1] - rcenter),xunitv)
        yunitv =   vec/norm(vec)  
        angles = []
        for i, vec in enumerate(facet):
            vx = dot(vec-rcenter,xunitv); vy = dot(vec-rcenter,yunitv)
            angle = arctan2(vy,vx)
            if angle < 0: angle = 2*pi + angle
            angles.append(angle)
        return [label for (angle,label) in sorted(zip(angles,labels))]

def flatVecsList(vecsList):
    '''Returns a flat list of vectors, from a list of lists that might have duplicates, 
    as in the case of facets on a polyhedron'''
    allpoints = []
    for isub,sub in enumerate(vecsList):
        for vec in sub:
            addVec(vec,allpoints)
    return allpoints

def plane3pts(points):
    '''From the first 3 points of a list of points, 
    returns a normal and closest distance to origin of the plane
    The closest distance is d = dot(r0,u).  The normal's direction 
    (degenerate vs multiplication by -1)
    is chosen to be that which points away from the side the origin is on.  If the 
    plane goes through the origin, then all the r's are in the plane, and no choice is made '''
    r0 = points[0]; r1 = points[1]; r2 = points[2]
    vec = cross(r1-r0,r2-r0)
    u = vec/norm(vec)
    if dot(u,r0) < 0: u = -u
    return u,dot(r0,u) 

def nextpos(ip,direc,nfacet):
    ip += direc
    if ip < 0: ip  += nfacet
    elif ip > nfacet: ip  -= nfacet
    return ip
    
class meshConstruct(): 
    '''Compact mesh reduced by point symmetry operations'''
    from comMethods import readfile,writefile,trimSmall,areEqual,directFromCart,cartFromDirect
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
    from conjGradMin2 import (fmin_cg,minimize_cg,line_search_wolfe1,scalar_search_wolfe1)

    def __init__(self):
        self.strLen = 5
        self.indIn = []
        self.indOut = []
        self.icut = -1
        
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
        self.eps = self.ravg/1000
        self.getVorCell() ############################   
        self.getIBZ() 
        self.meshCubic('fcc')
#         self.triFaces()
#         self.meshCubic('bcc')

#         self.meshCubic('cub')   
        
        sys.exit('stop')        
        print 'Typical spacing:', self.ravg
        self.initPoints() 
        self.initTets()
        self.movetoVoronoi()
        self.addBorders()
        self.delauPoints()
        self.fillTets()
        self.relax()
        self.plotPoints(self.points,self.npoints,'vec')
        sys.exit('stop')     
        return meshvecs, Nmesh, lattype, pfB, pf, status

#     def triFaces(self):
#         '''cell surface meshing'''
#         self.surfMeshPoints = []
#         for ifac,facet in enumerate(self.fpoints):
#             #Find 2-d representation for the points in this plane
#             inPlane = []
#             rcenter = sum(facet)/len(facet)
#             xunitv =  (facet[0] - rcenter)/norm(facet[0] - rcenter)
#     #         yunitv = cross(xunitv, rcenter)/norm(cross(xunitv, rcenter))    
#             vec = (facet[1] - rcenter) - dot((facet[1] - rcenter),xunitv)
#             yunitv =   vec/norm(vec)  
#             for i, vec in enumerate(facet):
#                 inPlane.append([dot(vec-rcenter,xunitv),dot(vec-rcenter,yunitv)])  
#             tri = delaunay(array(inPlane))
#             self.surfMeshPoints.append(delaunay())  
#         print len(self.surfMeshPoints),self.surfMeshPoints
#         self.allMathPrint()
    
    def meshCubic(self,type):
        '''Add a cubic mesh to the volume. If any 2 or 3 of the facet planes are 
        orthogonal, align the cubic mesh with them.       
        Weighting of each point:
        Volume of IBZ:  http://scipy.github.io/devdocs/generated/scipy.spatial.ConvexHull.html

        An "unprojected" point on the surface of IBZ: apply all symmetry operators.  
        Count how many new points it generates. 

        A projected point.  apply all symmetry operators.  
        Count how many new points it generates.  Its weight should then be 
        devalued vs an unprojected point on the surface by how much volume belongs to it.  
        We make a sphere for each point, which for normal packing is a given for a cell.  
        We derate the weight by how much the projected point's sphere overlaps other spheres.
        The volume of a sphere cut by a place a distance d<R from the center is 
        V = 1/3 Pi d(R-d)^2 (MathWorld).  d = ro-dot(r,u)

        More accurate: each point inside and across the border has its own mesh (e.g. bcc or fcc) voronoi cell. 
        If we keep (calculate) kpoints that have any of their facets inside the IBZ planes, 
        we weight them by their volume that is inside the 1BZ, even though the center might be outside.
        This can be done with machine precision.

        In practice, if a mesh point is within a distance of the packing sphere 
        radius from a plane, we use that plane to slice the VC of the point.   The distance is ro-dot(r,u).  Then use spatial.ConvexHull
        to get the volume'''
        a = 1.0
        cubicLVs = identity(3)
        if type == 'fcc':
            sites = [array([0, 0 , 0]), array([0,a/2,a/2]), array([a/2,0,a/2]), array([a/2,a/2,0])]
            pf = 0.74
        elif type == 'bcc':
            sites = [array([0, 0 , 0]), array([a/2,a/2,a/2])]
            pf = 0.68
        elif type == 'cub':
            pf = 0.52
            sites = [array([0, 0 , 0])]
        else: 
            sys.exit('stop. Type error in meshCubich')
        #test vacet points for orthogonality
        points = flatVecsList(self.fpoints)
        bodyCenter = sum(points)
        rs = []
        pairs = []
        triples = []
        for i in range(len(points)):
            if areEqual(norm(points[i]),0.0): break
            rs.append(norm(points[i]))
            for j in range(i,len(points)):
                if areEqual(norm(points[j]),0.0): break
                if areEqual(dot(points[i],points[j]),0.0):
                    pairs.append([points[i],points[j]])

        for ip,pair in enumerate(pairs):
            for i in range(len(points)):
                if areEqual(norm(points[i]),0.0): break
                if areEqual(dot(pair[0],points[i]),0.0) and areEqual(dot(pair[1],points[i]),0.0):
                    triple = deepcopy(pair)
                    triple.append(points[i])
                    triples.append(triple)
                    break
        
        #Define basis vectors for cubic lattice:
        Lsum= [] #length of vectors in pair or triplet
        if len(triple)>0:
            print 'At least one triplet of orthogonal point vectors found:',triples[0]
            if len(triples)>1: #find the one with most total vector length
                sums = zeros(len(triples))
                for it, triple in enumerate(triples):
                    for i in range(3):
                        sums[it] += norm(triple[i])
                triples = [triple for (sum1,triple) in zip(sums,triples)] #sorted by lowest sum    
            for i in range(3):
                vec = triples[-1][i]
                cubicLVs[:,i] = vec/norm(vec)
        elif len(pairs)>0:
            print 'At least one pair of orthogonal point vectors found:', pairs[0]
            if len(pairs)>1:
                sums = zeros(len(pairs))
                for ip, pair in enumerate(pairs):
                    for i in range(2):
                        sums[it] += norm(triple[i])
                pairs = [pair for (sum1,pair) in zip(sums,pairs)] #sorted by lowest sum    
            for i in range(2):        
                vec = pair[-1][i]
                cubicLVs[:,i] = vec/norm(vec)
            cubicLVs[:,2] = cross(cubicLVs[:,0],cubicLVs[:,1])
        else:
            print 'no orthogonal point vectors found.'
        #scale
        volSC= det(self.B)
        volKcubConv = volSC/self.nTarget*len(sites)
        aKcubConv = volKcubConv**(1/3.0)
        #Find the planar boundaries
        self.boundaries =[[],[]] #planes [u's],[ro's]
        for ifac, facet in enumerate(self.fpoints):
            facCenter = sum(facet)
            u,ro = plane3pts(facet[:3])
            if ro == 0: #plane goes through origin...choose a normal that points "outside" the cell
                if dot(facCenter - bodyCenter, u )<0:
                    u = -u
            self.boundaries[0].append(u);self.boundaries[1].append(ro)
        #Find the extremes in each cubLV direction:
        intMaxs = [] #factors of aKcubConv
        intMins = []
        for i in range(3):
            projs = []
            for point in points:
                projs.append(dot(cubicLVs[:,i],point))
            intMaxs.append(int(ceil(max(projs)/aKcubConv)))
            intMins.append(int(floor(min(projs)/aKcubConv)))       
        #Create the cubic mesh inside the irreducible BZ
        cubicLVs = cubicLVs * aKcubConv
        sites = [site * aKcubConv for site in sites]
        self.cubPoints = []
        self.cubWeights = []
#         offset = array([0.5,0.5,0.5])*aKcubConv
        for i in range(intMins[0],intMaxs[0]):
            for j in range(intMins[1],intMaxs[1]):
                for k in range(intMins[2],intMaxs[2]):
                    lvec = i*cubicLVs[:,0]+j*cubicLVs[:,1]+k*cubicLVs[:,2]
                    for site in sites:
                        kpoint = lvec + site
                        if not self.isOutside(kpoint):
                            self.cubPoints.append(kpoint)
                            self.cubWeights.append(self.IBZvolCut) #weight of a general point
        # We need to make our cubic mesh extend 
        # A mesh *point* is either outside, inside, on surface, on segment beween verticies, or is a vertex.
        # A point's volume may be inside, partially inside.   
        # If a point's sphere is cut by one plane, use the sphere's inside volume ratio to weight the 
        # point compared to a wholly inside point. This is fast.
        # If a point's sphere is cut by two or three planes (only a relatively few points),
        # use the VC cutting routines to find the volume inside.  
        # This automatically gives us all the symmetry weights!
                    
#         print 'cubPoints',len(self.cubPoints), self.cubPoints  
        self.facetsCubMathPrint(); print 'Show[p,q]'
        return     
        
    def magGroup(self,arr, igroup):
        '''returns the ith group, numbered from 1, (starting index and length of group) of vectors in a structured array, 
        which must be sorted by magnitude (in either direction), in field 'mag' '''

        newMags = [0] # shows indices of new group starts
        ngroups = 1
        for i in range(1,len(arr)):
            if abs(arr[i]['mag']-arr[i-1]['mag']) > self.eps:
                newMags.append(i)
                ngroups += 1
                if ngroups > igroup:
                    return newMags[ngroups-2],newMags[ngroups-1]-newMags[ngroups-2]
        return 0,len(arr)       
        
    def checkInside(self,grp):
        newPlanes = False
        newboundaries = self.boundaries
        for ig  in range(len(grp)):
            gvec = self.braggVecs[grp[ig]]['vec']
            if self.isInside(gvec):
                gnorm = norm(gvec)
                newboundaries[0].append(array(gvec/gnorm)); newboundaries[1].append(gnorm)
                newPlanes = True         
        self.boundaries = newboundaries
        return newPlanes
    
    def isInside(self,vec):
        '''Inside means on opposite side of the plane vs its normal vector'''
        inside = zeros(len(self.boundaries[0]),dtype = bool)
        for iplane, uvec in enumerate(self.boundaries[0]):
            if dot(vec,uvec) < self.boundaries[1][iplane] - self.eps: #point is inside this plane
                inside[iplane] = True
        return all(inside)

    def isOutside(self,vec):
        #print 'intersection',vec
        for iplane, uvec in enumerate(self.boundaries[0]): 
            pvec = uvec*self.boundaries[1][iplane]           
            if dot(vec,uvec) > self.boundaries[1][iplane] + self.eps: #point is outside this plane
#                 print '\tboundary', iplane, uvec, norm(pvec)
#                 print '\tboundary component', dot(vec,pvec), dot(pvec,pvec)
#                 print
                return True
#             else:
                #print 'On boundary check',iplane,dot(vec,pvec) - dot(pvec,pvec),pvec
        #print
        return False
              
    def onPlane(self,vec,planevec):
        return abs(dot(vec,planevec) - dot(planevec,planevec)) < self.eps #point is inside this plane
        
    def getVorCell(self):
        '''Boundaries and vertices of Voronoi cell'''
        self.getBraggVecs()
        igroup = 1
        checkNext = True
        gstart,ng = self.magGroup(self.braggVecs,1) # group of smallest bragg plane vectors
        self.boundaries =[[],[]] #planes [u's],[ro's]
        for i in range(ng):
            vec = self.braggVecs[i]['vec']; mag = norm(vec)
            self.boundaries[0].append(vec/mag); self.boundaries[1].append(mag)
        #print 'Smallest bragg vectors  boundaries', self.boundaries
        while checkNext:
            igroup += 1
            gstart,ng = self.magGroup(self.braggVecs,igroup)
            nextGroup = range(gstart,gstart+ng)
            checkNext = self.checkInside(nextGroup)
        self.getInterscPoints(self.boundaries)
        #write plane equations for Mathematica:
#         for iplane, uvec in enumerate(self.boundaries[0]):
#             pvec = uvec*self.boundaries[1][iplane]
#             print '{}x+{}y+{}z<={}&&'.format(pvec[0],pvec[1],pvec[2],dot(pvec,pvec)) #write plane equations for mathematica
        #print 'intersections'
#         for i in range(self.nIntersc):
            #print i, self.interscPoints[i]['mag'], self.interscPoints[i]['vec']
        self.getFacetsPoints()
        #print 'facet points'
#         for i in range(len(self.facetsPoints)):
#             print i, self.facetsPoints[i]['mag'], self.facetsPoints[i]['vec']
        self.arrangeFacets()
        return

    def arrangeFacets(self):        
        #arrange facets in each plane according to their angular order in the plane
        self.facets = [[]]*len(self.boundaries[0])
        for iplane, uvec in enumerate(self.boundaries[0]):
            pvec = uvec*self.boundaries[1][iplane]
            facetvecs = []
            facetlabels = []
            for i, fvec in  enumerate(self.facetsPoints['vec']):
                if self.onPlane(fvec,pvec):
                    facetlabels.append(i)
                    facetvecs.append(fvec)          
            self.facets[iplane] = orderAngle(facetvecs,facetlabels)
            #print iplane, 'facets', self.facets[iplane]
        
    def getFacetsPoints(self):
        OKpoints = []
        for i in range(len(self.interscPoints)):
            #print i,
            if not self.isOutside(self.interscPoints[i]['vec']):
                OKpoints.append(i)
        self.facetsPoints = zeros(len(OKpoints),dtype = [('vec', '3float'),('mag', 'float')])
        for iOK,ipoint in enumerate(OKpoints):
            vec = self.interscPoints[ipoint]['vec']
            self.facetsPoints[iOK]['vec'] = vec
            mag = norm(vec)
            self.facetsPoints[iOK]['mag'] = mag

    def getInterscPoints(self,planes):
        '''intersection points of planes, taken 3 at a time, where the planes
        have unit vector u (list 0) and are displaced a distance ro (list 1 in planes)'''
        vecs = [planes[0][i]*planes[1][i] for i in range(len(planes[0]))]
        combs = list(combinations(vecs,3))
        unique = []
        uniqStrs = []
        for c in combs:
            interscP = threePlaneIntersect(c)
            if not interscP is None:
                vecStr = str(interscP[0])[:self.strLen]+str(interscP[1])[:self.strLen]+str(interscP[2])[:self.strLen]
                if vecStr not in uniqStrs:
                    uniqStrs.append(vecStr)
                    unique.append(interscP)
        self.nIntersc = len(unique)
        self.interscPoints = zeros(self.nIntersc,dtype = [('vec', '3float'),('mag', 'float')])
        #print 'unique:',self.nIntersc
        for ipoint, vec in enumerate(unique):
            self.interscPoints[ipoint]['vec'] = vec
            self.interscPoints[ipoint]['mag'] = norm(vec)
#             print vec
        self.interscPoints.sort(order = 'mag')    
                        
    def getBraggVecs(self):
        '''The Bragg vectors are halfway from the origin to a lattice point.
        The planes normal to some of these will be bounding planes of the Voronoi cell '''
        self.braggVecs = zeros(124,dtype = [('vec', '3float'),('mag', 'float'),('dep', 'S15')])        
        ipoint = 0
        for i in range(-2,3):
            for j in range(-2,3):
                for k in range(-2,3):
                    if not (i==0 and j==0 and k==0):
                        vec = trimSmall(0.5*dot(self.B,array([i,j,k])))
                        self.braggVecs[ipoint]['vec'] = vec
                        self.braggVecs[ipoint]['dep'] = '{},{},{}'.format(i,j,k)
                        self.braggVecs[ipoint]['mag'] = norm(vec)
                        ipoint+=1
        self.braggVecs.sort(order = 'mag')
    
    def choose111(self,uvec):
        if dot(uvec,array([1,1,1]))>= 0:
            return trimSmall(real(uvec))
        else:
            return trimSmall(-real(uvec)) 
    
    def IntsPlLinSeg(self,u,r1,r2):
        '''Intersection between a plane through the origin and a line.
        A plane through the origin is given by dot(r,u) = 0.  
        A line segment between r1 and r2 is given by vectors r = r1+t(r2-r1), for t in (0,1) 
        Combining these:  r1u = dot(r1,u).  r2u = dot(r2,u).  Then t = r1u/(r1u-r2u).
        So if t is in (0,1), then we have an intersection'''
        den = dot(r1-r2,u)
        if areEqual(den,0.0):
            intersect = False
            return intersect, array([1e6,1e6,1e6])
        else:
            t = dot(r1,u)/den
            if 0 + self.eps < t < 1.0-self.eps:
                rinters =  trimSmall(r1 + t*(r2-r1))
                intersect = True
            else:
                 intersect = False
                 rinters =   array([1e6,1e6,1e6])
        return intersect, rinters        
    
    def cutCell(self,u):
        '''Cell is cut about the plane given by normal u.  Facets that intersect
        the plane are cut, and only the portion on one side is kept.  The intersection points
        between the plane and the facet segments are new facet points.  If a facet
        point lies on the plane, it stays in the facet. '''
        #print '\nu',u
#         if abs(norm(u-array([ 0.70710678,  0.,         -0.70710678]))) <1e-4:
#             print
        allRemoved = [] #points that are cut out
        bordersFacet = [] #new facet from the new edges of cut facets
#         ftemp = [[]]*len(self.fpoints) #this will contain only points, not labels
        ftemp = deepcopy(self.fpoints) #this will contain only points, not labels
#         if allclose(u,array([ 0.    ,      0.70710678 ,-0.70710678])):
#             print
        
        for ifac, facet in enumerate(self.fpoints):
            marker = ''
            #print 'facet',ifac, 'len',len(facet), facet
            bounds = []
            rbounds = []
            allLs = range(len(facet))
            keepLs = []
            for ip,point1 in enumerate(facet):
#                     marker = "*"                  
                if areEqual(dot(u,point1),0.0): #then this point is on the cut plane
                    bounds.append(ip)
                    rbounds.append('') #placeholder
                else: #check line segment
                    jp = ip + 1
                    if jp == len(facet): jp = 0 
                    point2 = facet[jp]
                    [intersect, rinters] = self.IntsPlLinSeg(u,point1,point2)          
                    if intersect:
                        bounds.append(ip + 0.5)  #to show it intersects, but not necessarily at a midpoint. 
                        rbounds.append(rinters)
           
            if 1 < len(bounds) < len(facet):
                if (isinteger(bounds[0]) and not (bounds[1]-bounds[0] ==1 or  
                        bounds[1]==len(facet)-1 and bounds[0] == 0) )\
                or\
                   (not isinteger(bounds[0]) or not isinteger(bounds[1])) : #if adjacent facet points, then one edge of a facet is on the plane...don't cut
                    #we check the u-projection of the point just beyond the first 
                    # boundary to see which part of the facet we keep...
                    # the half that u points away from
                    #print 'bounds',bounds
                    ucheck = dot(u,facet[int(ceil(bounds[0] + self.eps))])
                    direc = -int(sign(ucheck)) #+1 or -1 
                    newfacet = []
                    if isinteger(bounds[0]):
                        newfacet.append(facet[bounds[0]])
                        keepLs.append(bounds[0])
                        bordersFacet = addVec(facet[bounds[0]],bordersFacet)
                        #print 'append',bounds[0]
                        ip = nextpos(bounds[0],direc,len(facet))
                    else:
                        newfacet.append(rbounds[0])
                        bordersFacet = addVec(rbounds[0],bordersFacet)
                        #print 'append rb0'
                        ip = int(rint(bounds[0]+0.5*direc))
    #                 print ip,abs(ip - bounds[1])
                    while abs(ip - bounds[1])>=0.5 and abs(ip - (bounds[1] - len(facet)))>=0.5:
                        newfacet.append(facet[ip])
                        keepLs.append(ip)
                        #print 'append',ip
                        if abs(ip - bounds[1])==0.5 or abs(ip - (bounds[1] - len(facet)))==0.5: #next to boundary...done interior point
                            break
                        ip = nextpos(ip,direc,len(facet))
                    if isinteger(bounds[1]):
                        newfacet.append(facet[bounds[1]])
                        bordersFacet = addVec(facet[bounds[1]],bordersFacet)
                        keepLs.append(bounds[1])
                        #print 'append',bounds[1]
                    else:
                        newfacet.append(rbounds[1])
                        bordersFacet = addVec(rbounds[1],bordersFacet)
                        #print 'append rb1'
                    if len(newfacet) != len(facet) or not allclose(newfacet,facet):
                        #print 'newfacet',marker, newfacet
                        ftemp[ifac] = newfacet
                        #mark removed points for deletion in other facets
                        removed = [facet[i] for i in allLs if i not in keepLs]
                        for point in removed:
                            addVec(point,allRemoved)      
        if len(allRemoved)>0:
            self.icut += 1
            ftemp2 = deepcopy(ftemp)
            for i2, facet in enumerate(ftemp):
                nextbreak = False
                for ip,point in enumerate(facet):
                    for rpoint in allRemoved: 
                        if allclose(point,rpoint):
                            ftemp2[i2] = []
                            nextbreak = True
                            break
                    if nextbreak:
                        break
            ftemp = []
            for i2, facet in enumerate(ftemp2):
                if len(facet)> 0:
                   ftemp.append(facet)
            #Add any points that are in the cut plane into bordersFacet.  Some are not in facets with cuts. 
            points = flatVecsList(self.fpoints)
            for i in range(len(points)):
                if areEqual(dot(u,points[i]),0):
                    addVec(points[i],bordersFacet)            
            #Order by angle in the facet
            if len(bordersFacet)> 0:
                labels = orderAngle(bordersFacet,range(len(bordersFacet)))
                ftemp.append([bordersFacet[i] for i in labels])
                self.fpoints = ftemp
                #print 'Cut', self.icut
                self.facetsMathPrint(); print 'Show[p]'              
        return
                                              
#                 if isinteger(bounds[0]):
#                     for ip in range(bounds[0],bounds[1]+1):
#                         newfacet.append(facet[ip])
#                 else: # the facet runs from rbounds[0]...rbounds[1], which are in between facet points

#                 else: # we go backwards through the facet
#                     if isinteger(bounds[0]):
#                         for ip in range(bounds[0],bounds[1]+1):
#                             newfacet.append(facet[ip])
#                     else: # the facet runs from rbounds[0]...rbounds[1], which are in between facet points
#                         newfacet.append(rbounds[0])
#                         for ip in range(bounds[0]+ 0.5,bounds[1]-0.51):
#                             newfacet.append(facet[ip])
#                         newfacet.append(rbounds[1])                    
                    
    def makesDups(self,op):
        '''Applies symmetry operator to all facet points. If any facet points are 
        moved on top of other facet points, then return true'''
        points = flatVecsList(self.fpoints)
        for i in range(len(points)):
            rpoint = dot(op,points[i])
            if allclose(points[i],rpoint):
                break #no rotation effect
            otherLabels = range(len(points))
            otherLabels.pop(i)
            for j in otherLabels:
                if allclose(rpoint,points[j]):
                    return True
        return False
           
    def getIBZ(self):
        '''
        Apply symmetry operators to a facet point O:, to get point P.
        as a convention, choose points and volumes 
        that have vertices with the highest sum of components
        Define the point O as the vertex with the highest sum of components (or one of them if degenerate).
        1. Inversion: slice half the volume away about any plane through the center.
           Choice: define plane O-InvO-P, where P is any neighboring point to O on a fact.
        2. Reflection: slice through reflection plane
            Det is -1.  Has one eigenvalue that is -1, and that is the plane normal.
            a.  the plane normal gives a unique plane centered at the origin.  Use this to cut cell. 
        3. Rotation: They have det=1, and only one eigenvalue that is 1, which
        corresponds to the rotation axis.   with two -1's
            a. If O is not on axis, axis + rO form a plane.  Cut cell in half.  Then use axis + RO as
            a cut plane.  Keep part that has O and P. 
            b. if O is on axis, choose a different point for O.  
        from the axis of rotation to some point O and to point RO
        4. Improper rotation (Sn):  inversion + rotation, equivalent to rotation and reflection.
           Det -1.   In this case also we have a plane normal (reflection plane.  
           The cuts for improper rots are the same as for proper rots, so we change improper ones
           into proper ones before applying. 
        
        All rotations with angle not 0 or  pi have complex eigenvalues.  So an improper rotation
        has complex eigenvalues and a determinate less than 1. 
        
        
           
        Cutting: 
        Rotation: a plane through the origin has dot(r,u) = 0, for the normal u.  
        dot(r,u) > 0 means the half space beyond the plane in the direction of u. 
        We can get u defined by a point and the rotation axis  by u = +-cross(u_axis,u_point)
        For the +-, choose the direction to keep as that closest to the (1,1,1) direction
        (sum over components is largest). 
        
        Reflection/improper rotation: uR is given by the eigenvalue of R with eigen vector -1. 
        Choose either uR or -uR, the one closest to the (1,1,1) direction.   
        
        Keep the part of the cell that is dot(r,u) < 0

        '''
        #copy array of boundary plane normal vectors, that we will add to later
        

        #choose starting facet point arbitrarily by which has the highest sum of components
#         sumComps = []
        
#         for vec in self.facetsPoints[:]['vec']:
#             sumComps.append(sum(vec))
        print '\n\nReducing Brillouin zone by symmetry'
#         print '\n\nUsing dummy symmetry ops'
#         for iop in range(self.nops):
#             op = self.symops[:,:,iop]
#             print op        
        
#         self.newPlanes = []
        self.fpoints = [[]]*len(self.facets) #this will contain only point, not labels
        for ifac,facet in enumerate(self.facets):
            temp = []
            for ip in facet:
               temp.append(self.facetsPoints[ip]['vec'])
            self.fpoints[ifac] = temp
        print'Voronoi cell plot'; self.facetsMathPrint();print 'Show[p]' 
        inversion = False
        for iop in range(self.nops):
            op = self.symops[:,:,iop]            
            if abs(trace(op))< 3: #skip E and inverse
                evals,evecs = eig(op)
                evecs = array([evec for evec in evecs])
                if areEqual(det(op),-1.0) and not allclose(imag(evecs),zeros((3,3))):
#                     print '\nChanged improper rotation to proper one.'
                    op = -op
                    evals,evecs = eig(op)
                    evecs = array([evec for evec in evecs])
                #print '\niop',iop;print op
#                 if areEqual(det(op),1.0)  : #rotation
#                     print 'Rotation',  real(evecs[:,where(areEqual(evals,1.0))[0][0]])
#                 else:
#                     print 'Reflection',real(evecs[:,where(areEqual(evals,-1.0))[0][0]])
                if self.makesDups(op): #does this operation cause current facet points to lie on top of other current facet points. 
                    if areEqual(det(op),1.0)  : #rotation
                        evec = evecs[:,where(areEqual(evals,1.0))[0][0]] #axis
                        ipt = 0
                        #choose a facet point that is close to the rotation axis to avoid unusual cutting planes
                        ds = []
                        labels = range(len(self.facetsPoints))
                        for vec in self.facetsPoints['vec']:
                            if areEqual(abs(dot(evec,vec)),norm(vec)): #axis and vec are parallel...don't want this one.
                                 ds.append(100)
                            else:
                                ds.append(norm(vec - evec*dot(evec,vec)))
                        labels =  [label for (d,label) in sorted(zip(ds,labels))]#sort by distance
                        pnto = self.facetsPoints[labels[0]]['vec']
                        #the plane to cut is in the plane of O and axis, but perpendiular to vector O.                   
                        #tvec = cross(pnto,cross(pnto,evec))
                        tvec = cross(evec,pnto)
                        u1 = self.choose111(tvec/norm(tvec))
                        self.cutCell(u1)
                        pntp = dot(op,pnto)
                        tvec = cross(evec,pntp)
                        u2 = self.choose111(tvec/norm(tvec))
                        #print'\n\t2nd half of rot'
                        self.cutCell(u2) 

                    else: # -1: reflection/improper rotatioin
                        if len(where(areEqual(evals,-1.0))) > 1: evals = -evals #improper rotation
                        evec = evecs[:,where(areEqual(evals,-1.0))[0][0]]
                        u1 = self.choose111(evec) 
                        self.cutCell(u1) 
#                 else:
                    #print '\nOp {} yields no duplicates\n'.format(iop)  
            elif areEqual(det(op),-1.0):
                inversion = True
        if inversion: #apply last of all
            if self.makesDups(array([[-1.,  0.,  0.], [ 0., -1.,  0.], [ 0.,  0., -1.]])):
                #can cut along any plane
                self.cutCell(array([1.0,0.0,0.0]))
        #print 'facet points'
#         for ifac, facet in enumerate(self.fpoints):
#             print ifac, facet
        self.facetsMathPrint();print 'Show[p]'
        #get volume
        points = flatVecsList(self.fpoints)
        hull = convexH(points)
        self.IBZvolCut = det(self.B)/hull.volume
        print 'Vol BZ / Vol IBZ', self.IBZvolCut
        return

    def facetsMathPrint(self):
        ''' Mathematica'''
        print 'p = Graphics3D[{Red, Thick,{',
        for ifac, facet in enumerate(self.fpoints):
            print 'Line[{',
            for point in facet:
                print '{'+'{},{},{}'.format(point[0],point[1],point[2])+'},',
            print '{'+'{},{},{}'.format(facet[0][0],facet[0][1],facet[0][2])+'}', #back to first
            print '}]',
            if ifac < len(self.fpoints)-1:
                print ',',
        print '}}, Axes -> True,AxesLabel -> {"x", "y", "z"}]',    
        return    
    
    def facetsCubMathPrint(self):
        self.facetsMathPrint()
        print 'q=ListPointPlot3D[{',
        for ipoint,point in enumerate(self.cubPoints):
            print '{' + '{},{},{}'.format(point[0],point[1],point[2])+ '}',
            if ipoint < len(self.cubPoints) -1:
                print ',',
        print '},PlotStyle -> PointSize[0.03]];',
    
    def allMathPrint(self):
        ''' Mathematica'''
        self.facetsCubMathPrint()
        print 'r=ListPointPlot3D[{',
        for ipoint,point in enumerate(self.cubPoints):
            print '{' + '{},{},{}'.format(point[0],point[1],point[2])+ '}',
            if ipoint < len(self.surfMeshPoints) -1:
                print ',',
        print '},PlotStyle -> {Green,PointSize[0.03]}];',                
        print 'Show[p,q,r]'          
        return  
