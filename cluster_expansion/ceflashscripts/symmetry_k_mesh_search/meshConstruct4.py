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
from scipy.spatial import Delaunay as delaunay, Voronoi as sci_voronoi
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
    '''Makes superperiodic structures and cuts'''
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

    def meshCubic(self,type):
        '''Add a cubic mesh to the volume. If any 2 or 3 of the facet planes are 
        orthogonal, align the cubic mesh with them.  '''
        a = 1.0
        sitesBCC = [array([0, 0 , 0]), array([a/2,a/2,a/2])]
        sitesFCC = [array([0, 0 , 0]), array([0,a/2,a/2]), array([a/2,0,a/2]), array([a/2,a/2,0])]
        cubicLVs = identity(3)
        if type == 'fcc':
            sites = sitesFCC
            pf = 0.74
        elif type == 'bcc':
            sites = deepcopy(sitesBCC)
            pf = 0.68
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
        for i in range(intMins[0],intMaxs[0]):
            for j in range(intMins[1],intMaxs[1]):
                for k in range(intMins[2],intMaxs[2]):
                        lvec = i*cubicLVs[:,0]+j*cubicLVs[:,1]+k*cubicLVs[:,2]
                        for site in sites:
                            kpoint = lvec + site
                            if self.isInside(kpoint):
                                self.cubPoints.append(kpoint)
    
        print 'cubPoints',len(self.cubPoints), self.cubPoints  
        self.facetsMeshPrint(self.fpoints,self.cubPoints)
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
        print 'intersection',vec
        for iplane, uvec in enumerate(self.boundaries[0]): 
            pvec = uvec*self.boundaries[1][iplane]           
            if dot(vec,uvec) > self.boundaries[1][iplane] + self.eps: #point is outside this plane
#                 print '\tboundary', iplane, uvec, norm(pvec)
#                 print '\tboundary component', dot(vec,pvec), dot(pvec,pvec)
#                 print
                return True
            else:
                print 'On boundary check',iplane,dot(vec,pvec) - dot(pvec,pvec),pvec
        print
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
        print 'Smallest bragg vectors  boundaries', self.boundaries
        while checkNext:
            igroup += 1
            gstart,ng = self.magGroup(self.braggVecs,igroup)
            nextGroup = range(gstart,gstart+ng)
            checkNext = self.checkInside(nextGroup)
        self.getInterscPoints(self.boundaries)
        #write plane equations:
        for iplane, uvec in enumerate(self.boundaries[0]):
            pvec = uvec*self.boundaries[1][iplane]
            print '{}x+{}y+{}z<={}&&'.format(pvec[0],pvec[1],pvec[2],dot(pvec,pvec)) #write plane equations for mathematica
        print 'intersections'
        for i in range(self.nIntersc):
            print i, self.interscPoints[i]['mag'], self.interscPoints[i]['vec']
        self.getFacetsPoints()
        print 'facet points'
        for i in range(len(self.facetsPoints)):
            print i, self.facetsPoints[i]['mag'], self.facetsPoints[i]['vec']
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
            print iplane, 'facets', self.facets[iplane]
        
    def getFacetsPoints(self):
        OKpoints = []
        for i in range(len(self.interscPoints)):
            print i,
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
        print 'unique:',self.nIntersc
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
        print '\nu',u
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
            print 'facet',ifac, 'len',len(facet), facet
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
                    print 'bounds',bounds
                    ucheck = dot(u,facet[int(ceil(bounds[0] + self.eps))])
                    direc = -int(sign(ucheck)) #+1 or -1 
                    newfacet = []
                    if isinteger(bounds[0]):
                        newfacet.append(facet[bounds[0]])
                        keepLs.append(bounds[0])
                        bordersFacet = addVec(facet[bounds[0]],bordersFacet)
                        print 'append',bounds[0]
                        ip = nextpos(bounds[0],direc,len(facet))
                    else:
                        newfacet.append(rbounds[0])
                        bordersFacet = addVec(rbounds[0],bordersFacet)
                        print 'append rb0'
                        ip = int(rint(bounds[0]+0.5*direc))
    #                 print ip,abs(ip - bounds[1])
                    while abs(ip - bounds[1])>=0.5 and abs(ip - (bounds[1] - len(facet)))>=0.5:
                        newfacet.append(facet[ip])
                        keepLs.append(ip)
                        print 'append',ip
                        if abs(ip - bounds[1])==0.5 or abs(ip - (bounds[1] - len(facet)))==0.5: #next to boundary...done interior point
                            break
                        ip = nextpos(ip,direc,len(facet))
                    if isinteger(bounds[1]):
                        newfacet.append(facet[bounds[1]])
                        bordersFacet = addVec(facet[bounds[1]],bordersFacet)
                        keepLs.append(bounds[1])
                        print 'append',bounds[1]
                    else:
                        newfacet.append(rbounds[1])
                        bordersFacet = addVec(rbounds[1],bordersFacet)
                        print 'append rb1'
                    if len(newfacet) != len(facet) or not allclose(newfacet,facet):
                        print 'newfacet',marker, newfacet
                        ftemp[ifac] = newfacet
                        #mark removed points for deletion in other facets
                        removed = [facet[i] for i in allLs if i not in keepLs]
                        for point in removed:
                            addVec(point,allRemoved)      
        if len(allRemoved)>0:
            self.icut += 1
            if self.icut == 1:
                print
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
                print 'Cut', self.icut
                self.facetsMathPrint(self.fpoints)               
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
        print '\tThis routine does not reduce improper rotations'
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
        print'Voronoi cell plot'; self.facetsMathPrint(self.fpoints) 
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
                print '\niop',iop;print op
                if areEqual(det(op),1.0)  : #rotation
                    print 'Rotation',  real(evecs[:,where(areEqual(evals,1.0))[0][0]])
                else:
                    print 'Reflection',real(evecs[:,where(areEqual(evals,-1.0))[0][0]])
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
                        print'\n\t2nd half of rot'
                        self.cutCell(u2) 

                    else: # -1: reflection/improper rotatioin
                        if len(where(areEqual(evals,-1.0))) > 1: evals = -evals #improper rotation
                        evec = evecs[:,where(areEqual(evals,-1.0))[0][0]]
                        u1 = self.choose111(evec) 
                        self.cutCell(u1) 
                else:
                    print '\nOp {} yields no duplicates\n'.format(iop)  
            elif areEqual(det(op),-1.0):
                inversion = True
        if inversion: #apply last of all
            if self.makesDups(array([[-1.,  0.,  0.], [ 0., -1.,  0.], [ 0.,  0., -1.]])):
                #can cut along any plane
                self.cutCell(array([1.0,0.0,0.0]))
        print 'facet points'
        for ifac, facet in enumerate(self.fpoints):
            print ifac, facet
        self.facetsMathPrint(self.fpoints)

    def facetsMathPrint(self,farr):
        ''' Mathematica'''
        print 'Graphics3D[{Red, Thick,{',
        for ifac, facet in enumerate(farr):
            print 'Line[{',
            for point in facet:
                print '{'+'{},{},{}'.format(point[0],point[1],point[2])+'},',
            print '{'+'{},{},{}'.format(facet[0][0],facet[0][1],facet[0][2])+'}', #back to first
            print '}]',
            if ifac < len(farr)-1:
                print ',',
        print '}}, Axes -> True,AxesLabel -> {"x", "y", "z"}]'     
        return    
    
    def facetsMeshPrint(self,farr,parr):
        ''' Mathematica'''
        print 'p=Graphics3D[{Red, Thick,{',
        for ifac, facet in enumerate(farr):
            print 'Line[{',
            for point in facet:
                print '{'+'{},{},{}'.format(point[0],point[1],point[2])+'},',
            print '{'+'{},{},{}'.format(facet[0][0],facet[0][1],facet[0][2])+'}', #back to first
            print '}]',
            if ifac < len(farr)-1:
                print ',',
        print '}}, Axes -> True,AxesLabel -> {"x", "y", "z"}];',
        print 'q=ListPointPlot3D[{',
        for ipoint,point in enumerate(self.cubPoints):
            print '{' + '{},{},{}'.format(point[0],point[1],point[2])+ '}',
            if ipoint < len(self.cubPoints) -1:
                print ',',
        print '},PlotStyle -> PointSize[0.01]];', #end of ListPointPlot3D command 
        print 'Show[p,q]'          
        return  
  
    def movetoVoronoi(self):
        for i in range(self.npoints):
            self.points[i]['vec'] = intoVoronoi(self.points[i]['vec'],self.B)
    
    def volTet(self,vertPoints): 
        return abs(dot((vertPoints[0]-vertPoints[3]),cross((vertPoints[1]-vertPoints[3]),(vertPoints[2]-vertPoints[3])))/6.0)
    
    def delauPoints(self):
#         print 'len self.points', len(self.points[:self.npoints]['vec']) 
        tri = delaunay(self.points['vec'])
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
                self.tets[itet]['vol'] = abs(self.volTet([self.points[iv]['vec'] for iv in tet]))
        self.tets.sort(order='vol');self.tets = reverseStructured(self.tets)

    def addPntTet(self,itet,ipoint):
        tet = self.tets[itet]['tet'] 
        vp = array([self.points[iv]['vec'] for iv in tet])
        vd = array([self.points[iv]['dvec'] for iv in tet])
        newPoint = zeros(1,dtype = [('dep', 'S8'),('inCell', 'bool'),
            ('sponsor', 'int8'),('vec', '3float'),('dvec', '3float'),('force', '3float')])
        vec = sum(vp)/4.0 #center of mass
        newPoint[0]['vec'] = vec
#         print 'newvec',newPoint[0]['vec']
        dvec = sum(vd)/4.0 
        newPoint[0]['dvec'] = dvec
        newPoint[0]['dep'] = 'tet'.format(tet[0],tet[1],tet[2],tet[3]) 
        if min(dvec.flatten()) >= 0 and max(dvec.flatten()) < 1.0:
            newPoint[0]['inCell'] = True
            inCell = True
            self.nCellPoints += 1
            self.indIn.append(ipoint)
        else:
            newPoint[0]['inCell'] = False
            inCell = False
            self.indOut.append(ipoint)
        self.points = concatenate((self.points,newPoint),axis = 0)
        self.sympoints = [str(vec[0])[:self.strLen]+str(vec[1])[:self.strLen]+str(vec[2])[:self.strLen]]
        return vec,inCell 

    def addPntsSym(self,vec,isponsor,ipoint,type = None):
        '''Add points by point-symmetry.  Check that each is far enough away 
        from previous ones.  If too close, choose the point midway, make this the
        the independent point, and run through all the sym operators again.
        
        If the sponsor point is in the cell, we will move symmetry partners into cell.
        If not, we will not map them.  
        '''
        newPoints = zeros(self.nops -1,dtype = [('dep', 'S8'),('inCell', 'bool'),
                ('sponsor', 'int8'),('vec', '3float'),('dvec', '3float'),('force', '3float')])
        doneSym = False
        nNew = 0
        while not doneSym: #symmetry dependents
            vecList = [vec]
            ivec = 0
            for op in range(self.nops): 
                newvec = dot(self.symops[:,:,op],vec) #using voronoi cell. Sym ops don't take them out of cell. 
#                 if type == 'init':
#                     newvec = intoCell(dot(self.symops[:,:,op],vec),self.B)
#                 else: 
#                     newvec = dot(self.symops[:,:,op],vec)
                vecList.append(newvec)
                [tooClose,closevec] = self.tooCloseList(vecList,self.ravg/3.0)
                if not tooClose:
                    vecStr = str(newvec[0])[:self.strLen]+str(newvec[1])[:self.strLen]+str(newvec[2])[:self.strLen]
                    if vecStr not in self.sympoints:
                        self.sympoints.append(vecStr)
                        newPoints[ivec]['vec'] = newvec
                        dvecnew = self.directFromCart(self.B,newvec)    
                        newPoints[ivec]['dvec'] = dvecnew
                        newPoints[ivec]['dep'] = op
                        newPoints[ivec]['sponsor'] = isponsor
                        df = dvecnew.flatten()
                        inCell = (min(df) >= 0 and max(df) < 1 )
                        newPoints[ivec]['inCell'] = inCell
                        ivec += 1                          
                else:
                    #combine points, write over the last independent point, and start over on ops
                    print 'combine',newvec , closevec
                    vec = (newvec + closevec)/2.0
#                     if ipoint -1 in self.indIn: 
#                         self.indIn.remove(ipoint-1)
#                     else: 
#                         self.indOut.remove(ipoint-1)
                    self.points[ipoint-1]['vec'] = vec
                    dvec = self.directFromCart(self.B,vec)
                    self.points[ipoint-1]['dvec'] = dvec
                    self.points[ipoint-1]['inCell'] = inCell
                    self.sympoints = [str(vec[0])[:self.strLen]+str(vec[1])[:self.strLen]+str(vec[2])[:self.strLen]]
                    break                    
                    #keep dep the same: tet it came from...but this means that point may not be inside tet.
            doneSym = True
            for i in range(ivec):
                if newPoints[i]['inCell']: self.nCellPoints += 1                     
        self.points = concatenate((self.points,newPoints[:ivec]),axis = 0)
        return ivec
            


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
                [vec0,inCell] = self.addPntTet(0,ipoint)
                ipoint += 1
                isponsor = ipoint
                nNew = self.addPntsSym(vec0,isponsor,ipoint,inCell)
                ipoint += nNew       
                self.npoints = ipoint + 1 
                print'plotting {} points'.format(self.npoints)
                self.plotPoints(self.points,'vec',str(ipoint+1),nNew,self.tets[0:1])
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
                self.plotPoints(self.points,'vec',str(ipoint+1))   

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
                [vec0,inCell] = self.addPntTet(0,ipoint)
                ipoint += 1
                isponsor = ipoint
                nNew = self.addPntsSym(vec0,isponsor,ipoint,'init') #We are working in original parallelpiped
                ipoint += nNew   
#                 print'plotting {} points'.format(ipoint)
#                 self.plotPoints(self.points,'vec',str(ipoint),self.tets[0:1])
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
        self.plotPoints(self.points,'vec',str(self.nCellPoints),0,self.tets[:])
#         self.nCellPoints = ipoint + 1 - 8
        return       

    def initPoints(self):
        #self.subRandSym()
#         #use B cell vertices as the first 8 points
#         self.points = zeros(8,dtype = [('dep', 'S8'),('inCell', 'bool'),
#             ('sponsor', 'int8'),('vec', '3float'),('dvec', '3float'),('force', '3float')])
#         delta = 1e-6
#         a = 1-delta
#         z = 0+delta
#         for i,mult in enumerate([[z,z,z], [z,z,a],[z,a,z],[a,z,z],[z,a,a],[a,z,a],[a,a,z],[a,a,a]]):
#             dvec = array(mult)
#             self.points[i]['dvec'] = dvec
#             self.points[i]['vec'] = dot(self.B,dvec)
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
            ('sponsor', 'int8'),('vec', '3float'),('dvec', '3float'),('force', '3float')])
        
        self.points[0]['vec'] = array([0.0,0.0,0.0])
        ivec = 1
        for i in range(nvp):
            self.points[i]['vec'] = vorPoints[i]
            self.points[i]['inCell'] = True
            ivec += 1
        self.npoints = len(self.points)
        print'plotting {} points'.format(self.npoints)
        self.plotPoints(self.points,'vec',str(self.npoints))

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
                    for ivec in range(self.nCellPoints):
                        if not (i==0 and j==0 and k ==0):
                           transVec = dot(self.B,array([i,j,k]))
                           normTV = norm(transVec)
                           unitVec = transVec/normTV
                           newvec = self.points[ivec]['vec'] + transVec
                           if dot(newvec,unitVec) <= normTV/2.0 + self.rmax: #rmax gives the border thickness   
                               newPoint = zeros(1,dtype = [('dep', 'S8'),('inCell', 'bool'),
                                    ('sponsor', 'int8'),('vec', '3float'),('dvec', '3float'),('force', '3float')])
                               newPoint['dep'] = '{},{},{}'.format(i,j,k) # this point is dependent by translation
                               newPoint['sponsor'] = ivec #The sponsor may itself be dependent.  
                               newPoint['inCell'] = False                        
                               newPoint['vec'] = newvec
#                                self.points[ipoint]['dvec'] = newDirvec
                               self.points = concatenate((self.points,newPoint),axis = 0)
                               ipoint += 1
        self.npoints = ipoint
        self.points = delete(self.points,s_[self.npoints:],0)
        if len(self.points) != self.npoints: sys.exit('Stop.  Error in numbering points in expanded cell')
        print 'Init points plus border: {}'.format(self.npoints)
        self.plotPoints(self.points,'vec','expanded')
        return   
                    

        
#     def combineNear(self): 
#         '''Doesn't preserve symmetry, so not using it yet'''
#         self.cellPoints2 = zeros(self.nCellPoints+self.nops,dtype = [('dep', 'S8'),('sponsor', 'int8'),('vec', '3float'),('dvec', '3float'),('cell', '3int')])
# #         combine = []
#         for jvec in range(self.nCellPoints):
#             for kvec in range(jvec+1,self.nCellPoints):
#                 r = norm(self.cellPoints[jvec]['vec'] - self.cellPoints[kvec]['vec']) < self.rmin
#                 if 0.0 < norm(self.cellPoints[jvec]['vec'] - self.cellPoints[kvec]['vec']) < self.rmin:
#                     print jvec,kvec
#                     print 'r', norm(self.cellPoints[jvec]['vec'] - self.cellPoints[kvec]['vec'])
#                     if kvec in self.labels: self.labels.remove(kvec)
#                     if kvec in self.indIn: self.indIn.remove(kvec)
#                     if kvec in self.dep: self.dep.remove(kvec)
#                     self.cellPoints[jvec]['vec'] = (self.cellPoints[jvec]['vec'] + self.cellPoints[kvec]['vec'])/2.0     
#         ivec = 0
#         for label in self.labels:
#             self.cellPoints2[ivec] = self.cellPoints[label]
#             ivec += 1
#         self.nCellPoints = ivec
#         self.cellPoints[label] = self.cellPoints2[ivec]
#         print 'Points generated:',self.nCellPoints
#         return
                      
    def subRandSym(self):
        '''These points are in the VORONOI cell'''
#         self.labels = [] # list of active points labels
        self.cellPoints = zeros(self.nTarget+self.nops,dtype = [('dep', 'S8'),
            ('sponsor', 'int8'),('vec', '3float'),('dvec', '3float'),
            ('force', '3float')])
        self.directvec = zeros((self.nTarget,3))
        self.subrand = rand(3)
        afact = sqrt(2)
        bfact = sqrt(3)
        cfact = sqrt(5)
        self.indIn = []
        self.indOut = []
        ivec = -1
        edgeFactor = 10
        aEdge = self.rmax/norm(self.B[:,0])/edgeFactor; bEdge = self.rmax/norm(self.B[:,1])/edgeFactor; cEdge = self.rmax/norm(self.B[:,2])/edgeFactor        
        edges = [aEdge,bEdge,cEdge]
#         while ivec < self.nTarget: #do it once for now.
        while ivec < 0:
            itemp = ivec
            ivec += 1
            self.subrand = mod(self.subrand + array([afact,bfact,cfact]), array([1,1,1]))
            dvec = self.subrand
            vec = dot(self.B,dvec)
            vec = intoVoronoi(vec,self.B)#now in Voronoi cell
            print 'ind point', ivec,vec
            if ivec > 0 and self.tooCloseInit(vec,dvec,ivec,edges):
                print 'Independent point {} was too close. Trying again'.format(ivec)
                ivec = itemp
                continue #back to while
            self.cellPoints[ivec]['dep'] = 'I' #independent
            self.cellPoints[ivec]['sponsor'] = ivec                              
            self.cellPoints[ivec]['vec'] = vec 
            self.cellPoints[ivec]['dvec'] = dvec #dvec will always be the direct coordinates in the ORIGINAL cell
#             self.labels.append(ivec)
            self.indIn.append(ivec)
            iInd = ivec
            #now all symmetry operations will keep the point in the BZ 
            sympoints = [str(vec[0])[:self.strLen]+str(vec[1])[:self.strLen]+str(vec[2])[:self.strLen]]
            for op in range(self.nops):
                newvec = dot(self.symops[:,:,op],vec)
                vecStr = str(newvec[0])[:self.strLen]+str(newvec[1])[:self.strLen]+str(newvec[2])[:self.strLen]
                if vecStr not in sympoints:
                    dvecnew = self.directFromCart(self.B,newvec) 
                    if self.tooCloseInit(newvec,dvecnew,ivec,edges):
                        print 'Dependent point {} was too close. Trying again'.format(ivec)
                        ivec = itemp
                        self.indIn.remove(iInd)
                        break
                    sympoints.append(vecStr)
                    ivec += 1
#                     self.labels.append(ivec)
#                     self.dep.append(ivec)
                    self.cellPoints[ivec]['dep'] = '{}'.format(str(op))
                    self.cellPoints[ivec]['sponsor'] = iInd 
#                     newvec=intoVoronoi(newvec,self.B)                                
                    self.cellPoints[ivec]['dvec'] = dvecnew
#                     checkvec = intoVoronoi(newvec,self.B)
#                     if norm(checkvec-newvec)> 1e-6:
#                         print 'New Voronoi cell vec for',newvec,checkvec                   
                    self.cellPoints[ivec]['vec'] = newvec  
#                     print 'ivec',ivec, newvec 
#                     print '\tdvec',self.cellPoints[ivec]['dvec'] 
#                 if ivec == itemp: #some point is too close...try another independent point
#                     break                
        self.nCellPoints = ivec + 1
        self.cellPoints = delete(self.cellPoints,s_[self.nCellPoints:],0) #removing unused rows
        if len(self.cellPoints) != self.nCellPoints: sys.exit('Stop.  Error in numbering cellPoints')
        self.nIndIn = len(self.indIn) 
        self.plotPoints(self.cellPoints,'vec','voronoi')
#         self.plotPoints(self.cellPoints,'dvec')
        print 'Points in unit cell:',self.nCellPoints   
        print 'Independent points:',self.nIndIn                         
        return

    def tooCloseList(self,vecList,tol):
        '''Checks the last vec vs the others'''
        for partner in vecList[:len(vecList)-1]:
            r = norm(vecList[-1]-partner)
            if 1e-6 < r < tol:
                print 'Too close to another:',r
                return True, partner
        return False,[] 
    
    def tooCloseInit(self,vec,dvec,ivec,edges):
        for i in range(3):
            if abs(dvec[i]-0.5)> 0.5-edges[i]:
                print 'Too close to cell boundary'
                return True
        for testvec in self.cellPoints[:ivec]['vec']:
            if norm(vec-testvec) < self.ravg/5:
                print 'Too close to another:',norm(vec-testvec), testvec
                return True
        return False
    
    def plotPoints(self,pts,field = 'vec',tag = timestamp(),highlight = 0,tets = []):
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
        self.oldindVecs = [self.points[i]['vec'] for i in self.indIn]
        self.indInComps = array([self.points[i]['vec'] for i in self.indIn]).flatten()
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
                self.points[ipoint]['vec'] = intoCell(indVecs[place,:],self.B) 
#                 self.points[ipoint]['vec'] = indVecs[place,:]          
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
                    self.points[ipoint]['vec'] = intoCell(dot(self.symops[:,:,int(dep)],self.points[sponsor]['vec']),self.B)
                else: #translation
                    ms = array([int(i) for i in dep.split(',')])
                    svec = self.points[sponsor]['vec']
                    self.points[ipoint]['vec'] = svec + dot(self.B,ms)       
 
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
# #             print 'old',i,self.oldPoints[i]['vec']
# #             print 'new',i,self.points[i]['vec']
# #             print
#             move = norm(self.points[i]['vec']-self.oldPoints[i]['vec'])
#             if move > 1e-6 :
#                 print i,move
        
        ener = 0.0
        self.power = 2.0
#         scale = 1
        for ivec in range(self.nIndIn):
            force = zeros(3)
            rs = zeros(self.npoints,dtype = [('r', 'float'),('rij', '3float'),('force', '3float')])
                      
            nearest = -1
            rnearest = 100.0
            for jvec in range(self.npoints):
                rij = self.points[jvec]['vec'] - self.points[ivec]['vec']
                r = norm([rij])
                if r > 1e-4*self.rmin:
                    if r<rnearest:
                        nearest = jvec
                        rnearest = r
#                    
#                     epair = r**2 #attractive springs...dumb because those farthest away influence most
#                     force += rij #attractive springs)
                    
                    epair = (self.ravg/r)**self.power
                    forcepair = -self.power*epair/r * rij/r
                    force += forcepair
#                     print 'jvec',jvec,r,forcepair
                    ener += epair
                    rs[jvec]['r'] = r
                    rs[jvec]['rij'] = rij
                    rs[jvec]['force'] = forcepair
              
                    
            self.points[ivec]['force'] = trimSmall(force)
            
            if ivec in self.indIn:
                rs =sort(rs, order='r') 
#                 print ivec, rs[:20]['r'] 
#                 print rs[:20]['rij']
                print ivec, 'vec', self.points[ivec]['vec']
                print 'nerst',self.points[nearest]['vec'],nearest,rnearest 
                print 'vectr', trimSmall(self.points[nearest]['vec'] - self.points[ivec]['vec'])
                print 'force', ivec, self.points[ivec]['force']    
        ener = ener/self.nIndIn
        grad = -self.points[:self.nIndIn]['force'].flatten()/self.nIndIn 
        print 'energy:',self.count, ener
        print
#         if self.count == 19:
#             self.plotPoints(self.points,'vec')
        self.count += 1
        #now update all the dependent positions
      #  ! Note:  need to update the positions at each step!  Do we have to get inside
#         return ener#         
        return ener, grad
