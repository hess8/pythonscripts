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
    reverseStructured,isInVoronoi,areParallel, addVec, getSGpointGroup)

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
        
def orderAngle(facet,eps):
        '''get the angle that each vector is, in the plane of the facet.'''
        uvec = plane3pts(facet)[0] #normal of facet, not cut plane.  These are the same only for bordersfacet
        rcenter = sum(facet)/len(facet)
        xunitv =  (facet[0] - rcenter)/norm(facet[0] - rcenter)
        crss = cross(xunitv, uvec)
        yunitv = crss/norm(crss)
        angles = []
        for i, vec in enumerate(facet):
            vc = vec - rcenter
            vx = dot(vc,xunitv); vy = dot(vc,yunitv)
            angle = arctan2(vy,vx)
            if angle < 0-eps: angle += 2*pi
            angles.append(angle)
#         print 'angles',angles
        return [point for (angle,point) in sorted(zip(angles,facet),key = lambda x: x[0])] 

def flatVecsList(vecsList):
    '''Returns a flat list of vectors, from a list of lists that might have duplicates, 
    as in the case of facets on a polyhedron'''
    allpoints = []
    for isub,sub in enumerate(vecsList):
        for vec in sub:
            allpoints = addVec(vec,allpoints)
    return allpoints

def plane3pts(points):
    '''From the first 3 (noncolinear) points of a list of points, 
    returns a normal and closest distance to origin of the plane
    The closest distance is d = dot(r0,u).  The normal's direction 
    (degenerate vs multiplication by -1)
    is chosen to be that which points away from the side the origin is on.  If the 
    plane goes through the origin, then all the r's are in the plane, and no choice is made '''
    r0 = points[0]; r1 = points[1]; r2 = points[2]
    vec = cross(r1-r0,r2-r0)
    nv = norm(vec)
    if areEqual(nv,0):
        'Pause'
    u = vec/nv
    if dot(u,r0) < 0: u = -u
    return u,dot(r0,u)

def nextpos(ip,direc,nfacet):
    ip += direc
    if ip < 0: ip  += nfacet
    elif ip > nfacet: ip  -= nfacet
    return ip

def magGroup(arr, igroup,eps):
    '''returns the ith group, numbered from 1, (starting index and length of group) of vectors in a structured array, 
    which must be sorted beforehand by magnitude (in either direction), in field 'mag' '''
    newMags = [0] # shows indices of new group starts
    ngroups = 1
    for i in range(1,len(arr)):
        if abs(arr[i]['mag']-arr[i-1]['mag']) > eps:
            newMags.append(i)
            ngroups += 1
            if ngroups > igroup:
                return newMags[ngroups-2],newMags[ngroups-1]-newMags[ngroups-2]
    return 0,len(arr)   

def newBoundsifInside(braggVecs,grp,cell,eps):
    newPlanes = False
    newboundaries = deepcopy(cell.bounds)
    for ig  in range(len(grp)):
        gvec = braggVecs[grp[ig]]['vec']
        if isInside(gvec,cell.bounds,eps):
            gnorm = norm(gvec)
            newboundaries[0].append(array(gvec/gnorm)); newboundaries[1].append(gnorm)
            newPlanes = True         
    cell.bounds = newboundaries
    return newPlanes    

def getInterscPoints(planes):
    '''intersection points of planes, taken 3 at a time, where the planes
    have unit vector u (list 0) and are displaced a distance ro (list 1 in planes)'''
    vecs = [planes[0][i]*planes[1][i] for i in range(len(planes[0]))]
    combs = list(combinations(vecs,3))
    interscPoints = []
    interscPmags = []
    uniqStrs = []
    for c in combs:
        interscP = threePlaneIntersect(c)
        if not interscP is None:
            interscPoints = addVec(interscP,interscPoints)
    for ipoint, vec in enumerate(interscPoints):
        interscPmags.append(norm(vec))
    interscPoints = [point for (mag,point) in sorted(zip(interscPmags,interscPoints),key = lambda x: x[0])] #sorted by distance from origin 
    return interscPoints

def getFacetsPoints(interscPoints,cell,eps):
    cell.facets = [[]]*len(cell.bounds[0])
#     facetsPoints = []
#     for point in interscPoints:
#         if not isOutside(point,cell.bounds,eps):
#             facetsPoints.append(point)
    #arrange intersection points in each plane according to their angular order in the plane
    for iplane, uvec in enumerate(cell.bounds[0]):
        pvec = uvec*cell.bounds[1][iplane]
        facetvecs = []
        for i, vec in  enumerate(interscPoints):
            if onPlane(vec,pvec,eps) and not isOutside(vec,cell.bounds,eps):
                facetvecs = addVec(vec,facetvecs)          
        cell.facets[iplane] = orderAngle(facetvecs,eps)
    return cell

def isInside(vec,bounds,eps):
    '''Inside means on opposite side of the plane vs its normal vector'''
    inside = zeros(len(bounds[0]),dtype = bool)
    for iplane, uvec in enumerate(bounds[0]):
        if dot(vec,uvec) < bounds[1][iplane] - eps: #point is inside this plane
            inside[iplane] = True
    return all(inside)

def isOutside(vec,boundaries,eps):
    for iplane, uvec in enumerate(boundaries[0]): 
        pvec = uvec*boundaries[1][iplane]           
        if dot(vec,uvec) > boundaries[1][iplane] + eps: #point is outside this plane
            return True
    return False

def onPlane(vec,planevec,eps):
    return abs(dot(vec,planevec) - dot(planevec,planevec)) < eps #point is inside this plane
                    
def makesDups(op,cell):
    '''Applies symmetry operator to all facet points. If any facet points are 
    moved on top of other facet points, then return true'''
    points = cell.fpoints
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

def getVorCell(LVs,cell,eps):
    '''Boundaries and vertices of Voronoi cell'''
    braggVecs = getBraggVecs(LVs)
    igroup = 1
    checkNext = True
    gstart,ng = magGroup(braggVecs,1,eps) # group of smallest bragg plane vectors
    for i in range(ng):
        vec = braggVecs[i]['vec']; mag = norm(vec)
        cell.bounds[0].append(vec/mag); cell.bounds[1].append(mag)
    #print 'Smallest bragg vectors  boundaries', cell.bounds
    while checkNext:
        igroup += 1
        gstart,ng = magGroup(braggVecs,igroup,eps)
        nextGroup = range(gstart,gstart+ng)
        checkNext = newBoundsifInside(braggVecs,nextGroup,cell,eps)
    interscPoints = getInterscPoints(cell.bounds)
    #write plane equations for Mathematica:
#         for iplane, uvec in enumerate(cell.bounds[0]):
#             pvec = uvec*cell.bounds[1][iplane]
#             print '{}x+{}y+{}z<={}&&'.format(pvec[0],pvec[1],pvec[2],dot(pvec,pvec)) #write plane equations for mathematica
    #print 'intersections'
#         for i in range(self.nIntersc):
        #print i, interscPoints[i]['mag'], interscPoints[i]['vec']
    cell = getFacetsPoints(interscPoints,cell,eps)
    cell.fpoints = flatVecsList(cell.facets) 
    cell.center = sum(cell.fpoints)/len(cell.fpoints)
    cell.volume = convexH(cell.fpoints).volume
    return cell
    #print 'facet points'
#         for i in range(len(self.facetsPoints)):
#             print i, self.facetsPoints[i]['mag'], self.facetsPoints[i]['vec']
             
def getBraggVecs(LVs):
    '''The Bragg vectors are halfway from the origin to a lattice point.
    The planes normal to some of these will be bounding planes of the Voronoi cell '''
    braggVecs = zeros(124,dtype = [('vec', '3float'),('mag', 'float'),('dep', 'S15')])        
    ipoint = 0
    for i in range(-2,3):
        for j in range(-2,3):
            for k in range(-2,3):
                if not (i==0 and j==0 and k==0):
                    vec = trimSmall(0.5*dot(LVs,array([i,j,k])))
                    braggVecs[ipoint]['vec'] = vec
                    braggVecs[ipoint]['dep'] = '{},{},{}'.format(i,j,k)
                    braggVecs[ipoint]['mag'] = norm(vec)
                    ipoint+=1
    braggVecs.sort(order = 'mag')
    return braggVecs

def intsPlLinSeg(u,ro,r1,r2,eps):
    '''Intersection between a plane through the origin and a line.
    A plane through the origin is given by dot(r,u) = ro.  
    A line segment between r1 and r2 is given by vectors r = r1+t(r2-r1), for t in (0,1) 
    Combining these:  r1u = dot(r1,u).  r2u = dot(r2,u).  Then t = (ro-r1u)/(r2u-r1u).
    So if t is in (0,1), then we have an intersection'''
    r1u = dot(r1,u)
    r2u = dot(r2,u)
    if areEqual(r1u,r2u):
        return False, None
    else:
        t = (ro-r1u)/(r2u-r1u)
        if 0 + eps < t < 1.0-eps:
            return True, trimSmall(r1 + t*(r2-r1)) 
        else:
            return False, None    

def getBoundsFacets(cell,rpacking = None): 
    '''After cutting, we have new facets, and the facet planes are the boundaries'''
    cell.bounds = [[],[]] 
    
    for ifac, facet in enumerate(cell.facets):
        facCenter = sum(facet)/len(facet)
        u,ro = plane3pts(facet[:3])
        if ro == 0 and dot(facCenter - cell.center, u )<0: #plane goes through origin...choose a normal that points "outside" the BZ
            u = -u
        cell.bounds[0].append(u)
        cell.bounds[1].append(ro)
    if not rpacking is None:
        cell.expBounds = [[],[]]
        cell.expBounds[0] = cell.bounds[0]
        cell.expBounds[1] = [ro+sqrt(2)*rpacking for ro in cell.bounds[1]] #sqrt(2) for vertices of Voronoi cell vs radius of sphere.
    return cell

class cell():
    def __init__(self):
        self.bounds = [[],[]] #planes, written as normals and distances from origin [u's] , [ro's]
        self.expBounds = None #planes pushed outward by rpacking
        self.facets = None #points arranged in facets
        self.fpoints = None #points as a set (but a list)
        self.volume = None
        self.center = None #body center of cell
        self.mesh = None #centers of each voronoi cell, or kpoints
        self.weights = None
    
class meshConstruct(): 
    '''Compact mesh reduced by point symmetry operations'''
    from comMethods import readfile,writefile,trimSmall,areEqual,directFromCart,cartFromDirect
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
    def __init__(self):
        self.icut = -1
#         self.BZboundaries = [[],[]] #for BZ Voronoi cell or IBZ.  written as [u's] , [ro's]
#         self.MPboundaries = [[],[]] #for mesh point Voronoi cell
#         cell.facets = [] #facet points for BZ Voronoi cell or IBZ
#         self.MPfacets = [] #facet points for mesh point Voronoi cell
        
    def meshSym(self,A,B,totatoms,aTypes,postype,aPos,targetNmesh,path,method):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        self.B = B
        print 'B',B
        [self.symops,self.nops] = getGroup(self.B)
#         [self.symops,self.nops] = getGroup(A)
#         [self.symops,self.nops] = getGroup(transpose(A))
#         [self.symops,self.nops] = testFor(A,aTypes,aPos)
#         if postype[0].lower() == 'd': #direct
#             direct = True
#         else:
#             direct = False
#         [self.symops,self.nops] = getSGpointGroup(A,totatoms,aTypes,aPos,direct)
#         [self.symops,self.nops] = getSGpointGroup(transpose(A),aTypes,aPos)
        self.nTarget = targetNmesh
#         self.nTarget = 48
#         print 'testing with ntarget',self.nTarget
        self.path = path
        self.method = method
            #0,1 for now.
            #0: exact: use vertices of mesh voronoi cell that are closest/farthest 
            #         from the IBZ center origin to check if the point's volume is cut. 
            #         Cut the VC to determine the volume contribution      
            #1: approx 1. Use sphere around mesh point to test whether it is near a surface.  
            #         For a 1-plane cut, use the spherical section that is inside. 
            #         For 2 or 3 plane cut, we use the exact method.
        print 'Number of desired points:', targetNmesh
        print 'Symmetry operations:', self.nops
        vol = abs(det(B))
        self.ravg = (vol/targetNmesh)**(1/3.0)
        self.edgeFactor = 3.0
        self.rmin = 0.8*self.ravg #
        self.rmax = self.edgeFactor*self.ravg #cutoff for forces, neighbors in other cells. 
        eps = self.ravg/1000
        BZ = cell() #instance
        BZ = getVorCell(self.B,BZ,eps)
        self.facetsMathPrint(BZ,'p',True,'Red') 
        BZ = self.getIBZ(BZ,eps) #changes BZboundaries.         
        self.meshCubic(BZ,'fcc',eps)
#         self.triFaces()
#         self.meshCubic('bcc')

#         self.meshCubic('cub')   
        self.writeKpoints(BZ)       
        return
    
    def writeKpoints(self,cell):
        nk = len(cell.mesh)
        totw = sum(cell.weights)
        lines = []
        lines.append('Vornoi cell tiling of IBZ (Bret Hess, BYU). Total weights: {:12.8f} (vs 1.0 per general point without symmetry\n'.format(totw))
        lines.append('{}\n'.format(nk))
        lines.append('Cartesian\n')
        for ik,kpoint in enumerate(cell.mesh):
            lines.append('{:15.12f}  {:15.12f}  {:15.12f}  {:15.12f}\n'\
                         .format(kpoint[0],kpoint[1],kpoint[2],cell.weights[ik]))
        self.writefile(lines,'KPOINTS')         
   
    def meshCubic(self,BZ,type,eps):
        '''Add a cubic mesh to the interior, . If any 2 or 3 of the facet planes are 
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
        The volume of a sphere cut by a plane a distance d<R from the center is 
        V = 1/3 Pi (R-d)^2(2R-d)
        More accurate: each point inside and across the border has its own mesh (e.g. bcc or fcc) voronoi cell. 
        If we keep (calculate) kpoints that have any of their facets inside the IBZ planes, 
        we weight them by their volume that is inside the 1BZ, even though the center might be outside.
        This can be done with machine precision.

        In practice, if a mesh point is within a distance of the packing sphere 
        radius from a plane, we use that plane to slice the VC of the point.   The distance is ro-dot(r,u).  Then use spatial.ConvexHull
        to get the volume
        
        # We need to make our cubic mesh extend one layer beyond any IBZ borders
        # A mesh *point* is either outside, inside, on surface, on segment beween verticies, or is a vertex.
        # A point's volume may be inside, partially inside.   
        # If a point's sphere is cut by one plane, use the sphere's inside volume ratio to weight the 
        # point compared to a wholly inside point. This is fast.
        # If a point's sphere is cut by wtwo or three planes (only a relatively few points),
        # use the VC cutting routines to find the volume inside.
        # This automatically gives us all the symmetry weights!'''
        
        a = 1.0
        cubicLVs = identity(3)
        if type == 'fcc':
            sites = [array([0, 0 , 0]), array([0,a/2,a/2]), array([a/2,0,a/2]), array([a/2,a/2,0])]
            primLVs = array([[0,a/2,a/2], [a/2,0,a/2], [a/2,a/2,0]])
            self.rpacking = 1/2.0/sqrt(2)
            pf = 4*4/3.0*pi*self.rpacking**3  #0.74
#             fccFacets = 
        elif type == 'bcc':
            sites = [array([0, 0 , 0]), array([a/2,a/2,a/2])]
            primLVs = array([[-a/2,a/2,a/2], [a/2,-a/2,a/2], [a/2,a/2,-a/2]])
            self.rpacking = sqrt(3)/4
            pf = 2*4/3.0*pi*self.rpacking**3 #0.68
        elif type == 'cub':
            sites = [array([0, 0 , 0])]
            primLVs = cubicLVs
            self.rpacking = 0.5
            pf = 4/3.0*pi*self.rpacking**3 #0.52
        else:
            sys.exit('stop. Type error in meshCubich')
        #test facet points for orthogonality
        rs = []
        pairs = []
        triples = []
        points = BZ.fpoints
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
        if len(triples)>0:
            print 'At least one triplet of orthogonal point vectors found:',triples[0]
            if len(triples)>1: #find the one with most total vector length
                sums = zeros(len(triples))
                for it, triple in enumerate(triples):
                    for i in range(3):
                        sums[it] += norm(triple[i])
                triples = [triple for (sum1,triple) in sorted(zip(sums,triples),key = lambda x: x[0])] #sorted by lowest sum  
            for i in range(3):
                vec = triples[-1][i]
                cubicLVs[:,i] = vec/norm(vec)
        elif len(pairs)>0:
            print 'At least one pair of orthogonal point vectors found:', pairs[0]
            if len(pairs)>1:
                sums = zeros(len(pairs))
                for ip, pair in enumerate(pairs):
                    for i in range(2):
                        sums[ip] += norm(pairs[i])
                pairs = [pair for (sum1,pair) in sorted(zip(sums,pairs),key = lambda x: x[0])] #sorted by lowest sum    
            for i in range(2):        
                vec = pair[-1][i]
                cubicLVs[:,i] = vec/norm(vec)
            cubicLVs[:,2] = cross(cubicLVs[:,0],cubicLVs[:,1])
        else:
            print 'no orthogonal point vectors pairs found.'
        #scale
        volSC= det(self.B)
        volKcubConv = volSC/self.nTarget*len(sites)
        aKcubConv = volKcubConv**(1/3.0)
        primLVs = primLVs * aKcubConv
        cubicLVs = cubicLVs * aKcubConv
        sites = [site * aKcubConv for site in sites]
        self.rpacking = self.rpacking*aKcubConv
        self.Vsphere = 4/3.0*pi*self.rpacking**3
        print 'rpacking',self.rpacking
        BZ = getBoundsFacets(BZ,self.rpacking) #adds expanded cell
        #define mesh point voronoi cell
        MP = cell()
        MP = getVorCell(primLVs,MP,eps)
        print 'Mesh Point Voronoi cell',MP.facets; self.facetsMathPrint(MP,'p',True,'Red');print ';Show[p]\n'        
        #Find the extremes in each cubLV direction:
        intMaxs = [] #factors of aKcubConv
        intMins = []
        for i in range(3):
            projs = []
            for point in points:
                projs.append(dot(cubicLVs[:,i],point)/aKcubConv**2)
            intMaxs.append(int(ceil(max(projs)))+1)
            intMins.append(int(floor(min(projs)))-1)       
        #Create the cubic mesh inside the irreducible BZ
        BZ.mesh = []
        BZ.weights = []
        weightsInside = 0
        nInside = 0
        weightsOnePlane = 0
        nOnePlane = 0
        weightsCuts = 0
        nCut = 0
#         offset = array([0.5,0.5,0.5])*aKcubConv
        print 'bounds'
        for i in range(len(BZ.bounds[0])):
            print BZ.bounds[0][i],BZ.bounds[1][i]               
        ik = 0 
        #begin MP facets printing
        self.facetsMathPrint(BZ,'s','True','Red'); print ';', #draw supecell voronoi cell
        showCommand = 'Show[s,' 
        #end start of MP facets printing
        for i in range(intMins[0],intMaxs[0]):
            for j in range(intMins[1],intMaxs[1]):
                for k in range(intMins[2],intMaxs[2]):
                    lvec = i*cubicLVs[:,0]+j*cubicLVs[:,1]+k*cubicLVs[:,2]
                    for site in sites:
                        ik+=1
                        kpoint = lvec + site
                        ds = self.dToPlanes(kpoint,BZ.expBounds)
                        if self.isAllInside(ds,eps):
#                             print '\nkpoint',kpoint, ik, [i,j,k]; 
#                             print 'ds',ds
#                             print 'inside',ds,kpoint
                            BZ.mesh.append(kpoint)
                            BZ.weights.append(self.IBZvolCut)
                            weightsInside += self.IBZvolCut
                            nInside += 1
                        elif self.isInsideExpanded(ds,eps):
                            #change to d's from real cell, the borders of the BZ
                            dsBZ = [d + sqrt(2)*self.rpacking for d in ds]
#                             print '\n\nkpoint',kpoint, ik, [i,j,k]; sys.stdout.flush()
#                             print 'ds',dsBZ
                            
                            
#                             if self.method == 0:
#                                 testMP = self.prepCutMP(MP,kpoint)
#                                 if ik == 19:
#                                     'pause'
#                                 nearPlanes,nearDs = self.checkVertices(BZ,testMP,dsBZ) #do any vertices lie beyond the BZ boundaries
#                                 if len(nearPlanes) == 0: 
#                                     break
#                             elif self.method == 1:
#                                 near = where(abs(dsBZ) <= 2.0*self.rpacking - eps)
#                                 nearPlanes = near[0] #arrays
#                                 nearDs = dsBZ[near]
#                             if len(nearPlanes)>0:
#                                 if len(nearPlanes) == 1 and self.method == 1: #give weight proportional to sphere
#                                     BZ.mesh.append(kpoint)
# #                                     d2 = ds[nearPlanes[0]] + self.rpacking #the d's are negative here, and measured vs the expanded boundary
# #                                     print 'one-plane cut',d2, kpoint
#                                     cap = 1/3.0*pi*(self.rpacking-abs(dsBZ[0]))**2 *(2*self.rpacking+abs(d2))
# #                                     print 'cap',cap, 'vs sphere',self.Vsphere
#                                     if d2<=0: #most of the sphere is inside:
#                                         weight = self.IBZvolCut*(1-cap/self.Vsphere)
#                                     else:
#                                         weight = self.IBZvolCut*cap/self.Vsphere
#                                     BZ.weights.append(weight)
#                                     weightsOnePlane += weight
#                                     nOnePlane += 1
# #                                     print 'weight',weight
#                                 elif len(nearPlanes) <=4: #cut the mesh point Voronoi cell:
# #                                     if ik == 54:
# #                                         self.facetsMathPrint(cutMP,'p',True,'Red');print ';Show[p]\n'
                            
                            cutMP = self.prepCutMP(MP,kpoint)  #cut MP is displaced by kpoint from MP
#                             if ik == 56:
#                                 'pause'
#                                     self.facetsMathPrint(cutMP,'p',True,'Red'); print ';Show[p]\n' 
#                                         'pause'
                            for iplane, uvec in enumerate(BZ.bounds[0]):
#                                 uvec = BZ.bounds[0][BZlabel]
                                ro = BZ.bounds[1][iplane]                                      
#                                 print 'u',uvec,'ro',ro
                                cutMP = self.cutCell(uvec,ro,cutMP,eps) # we always keep the part that is "inside", opposite u
#                                 self.facetsMathPrint(cutMP,'p',True,'Red'); print ';Show[p]\n' 
                                if cutMP.volume == 0.0: #(outside BZ. happens in oblique corners of expanded cell)
#                                     print 'volume is zero: zero points'
                                    break
#                                         if ik == 105:
#                                         self.facetsMathPrint(cutMP,'p',True,'Red');print ';Show[p]\n'
                                if len(cutMP.facets)>=4:
                                    cutMP.volume = convexH(cutMP.fpoints).volume
                                else:
                                    cutMP.volume = 0.0
#                                     print 'volume is zero: fewer than 4 facets'
                                    break
                            if cutMP.volume > 0.0:
                                BZ.mesh.append(kpoint)
                                weight = self.IBZvolCut*cutMP.volume/MP.volume; BZ.weights.append(weight)
    #                                    print 'weight',weight
                                weightsCuts += weight
                                nCut += 1
                                #####MP facet printing loop line
                                print 'ik',ik
#                                 showCommand = self.cutMPCellMathPrint(BZ,cutMP,kpoint,ik,showCommand)
                                #####end MP facet printing loop entry
        if showCommand[-1] == ',': showCommand = showCommand[:-1]
        showCommand += ']' 
        print ';', 
        print showCommand 
        #end MP facets printing
        print 'Weights inside', nInside, weightsInside
        if nOnePlane>0: print 'Weights one plane',nOnePlane,weightsOnePlane, weightsOnePlane/float(nOnePlane)
        print 'Weights cuts', nCut,weightsCuts, 'Average per cut VC', weightsCuts/float(nCut)
        print 'Total volume in weights:', (nCut + nOnePlane + nInside), sum(BZ.weights)*MP.volume
        print 'BZ volume:', det(self.B)
        self.facetsMeshMathPrint(BZ); print ';Show[p,q]\n'
        self.facetsMeshVCMathPrint(BZ,MP)
        return

    def prepCutMP(self,MP,kpoint):
#         cutMP = deepcopy(MP)
        cutMP = cell()
        cutMP.facets = [[]]*len(MP.facets)
        for ifac, facet in enumerate(MP.facets):
            temp = []
            for point in facet:
                temp.append(point + kpoint)
            cutMP.facets[ifac] = temp
        cutMP.fpoints = []
        for ipoint, point in enumerate(MP.fpoints):
            cutMP.fpoints.append(point + kpoint)
        cutMP.center = kpoint
        return cutMP
        
    def cutMPCellMathPrint(self,BZ,cellcut,point,ipoint,showCommand):
        '''For showing cut facets inside IBZ, loop work'''
#         tCell = cell()
#         tCell.facets = [[]]*len(cellcut.facets)
#         for ifac,facet in enumerate(cellcut.facets):
#             temp = []
#             for facetPoint in facet:
#                 temp.append(facetPoint + point)
#             tCell.facets[ifac] = temp
#         self.facetsMathPrint(cellcut,'v{}'.format(ipoint),False,'Blue');print ';',
        ncolors = 100
#         self.facetsMathPrint(cellcut,'v{}'.format(ipoint),False,'discreteColors[{}][[{}]]'\
#                              .format(ncolors,mod(ipoint,ncolors)));print ';',
        self.facetsMathPrint(cellcut,'v{}'.format(ipoint),False,'RandomColor[]');print ';',
        showCommand += 'v{},'.format(ipoint)
        return showCommand
            
    def dToPlanes(self,vec,bounds):
        ds = []
        for iplane, uvec in enumerate(bounds[0]):   
            d = dot(vec,uvec) - bounds[1][iplane] 
            ds.append(d)
        return array(ds)
    
    def choose111(self,uvec,eps):
        if dot(uvec,array([1,1,1]))> 0 + eps:
            return -trimSmall(real(uvec))
        else:
            return trimSmall(real(uvec))                   

    def cutCell(self,u,ro,cell,eps):
        '''Cell is cut about an arbitrary plane given by normal u and plane distance from 
        origin ro.  Facets that intersect the plane are cut, 
        and only the portion on one side is kept.  The intersection points
        between the plane and the facet segments are new facet points.  If a facet
        point lies on the plane, it stays in the facet.'''
        allRemoved = [] #points that are cut out
        bordersFacet = [] #new facet from the points of cut facets      
#         print 'u',u,'ro',ro
        ftemp = deepcopy(cell.facets)       
        for ifac, facet in enumerate(cell.facets):
#             print 'facet',ifac,'len',len(facet),facet
#             print 'facet',ifac, 'len',len(facet), facet
            newfacet = []
            #first, determine if the plane intersects any facet points or 
            #line segments, by checking the change of sign of dot products
#             sign0 = sign(dot(u,facet[0]))
            signs = []
            for point in facet:
                pu = dot(u,point)
                if pu > ro + eps:
                    signs.append(1.0)
                elif areEqual(pu,0.0):
                    signs.append(0.0)
                else:
                    signs.append(-1.0)
            if -1 in signs and 1 in signs: #facet is cut
                for ip,pointi in enumerate(facet):
                    pu = dot(u,pointi)
                    if areEqual(pu,ro): # then this point is on the cut plane 
                        newfacet.append(pointi)
                        bordersFacet = addVec(pointi,bordersFacet)
                    elif pu < ro + eps: #inside the plane
                        newfacet.append(pointi)
                    else: #outside
                        allRemoved = addVec(pointi,allRemoved)
                    #check line segment between point1 and next point on facet
                    jp = ip + 1
                    if jp == len(facet): jp = 0 
                    pointj = facet[jp]
                    [intersect, rinters] = intsPlLinSeg(u,ro,pointi,pointj,eps) #         
                    if intersect:
                        newfacet.append(rinters)
                        bordersFacet = addVec(rinters,bordersFacet)
#                 print 'newFacet',newfacet
#                 print 'bordersFacet',bordersFacet
                if len(newfacet) >= 3:
                    ftemp[ifac] = orderAngle(newfacet,eps)
            else: #mark for removal all points that are outside of the plane
                for i, sgn in enumerate(signs):
                    if sgn == 1.0:
                        allRemoved = addVec(facet[i],allRemoved)
                    elif sgn == 0.0:
                        bordersFacet = addVec(facet[i],bordersFacet)         
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
            #Add any points that are in the cut plane into bordersFacet.  Some may not be in facets with cuts. 
            points = flatVecsList(cell.facets)
            for point in points:
                if areEqual(dot(u,point),ro):
                    bordersFacet = addVec(point,bordersFacet)            
            #Order by angle in the facet
            if len(bordersFacet)> 2:
                uvec = plane3pts(bordersFacet)[0]
                ftemp.append(orderAngle(bordersFacet,eps))
            cell.facets = ftemp
            cell.fpoints = flatVecsList(cell.facets) 
            if len(cell.fpoints)== 0:
                cell.volume = 0
            else:
                cell.center = sum(cell.fpoints)/len(cell.fpoints)
#             self.facetsMathPrint(cell,'p',True,'Red');print ';Show[p]\n'              
        return cell      

    def getIBZ(self,BZ,eps):
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
        
        print '\n\nReducing Brillouin zone by symmetry'
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
                if makesDups(op,BZ): #does this operation cause current facet points to move to other facet points
                    if areEqual(det(op),1.0)  : #rotation
                        evec = evecs[:,where(areEqual(evals,1.0))[0][0]] #axis
                        ipt = 0
                        #choose a facet point that is close to the rotation axis to avoid unusual cutting planes
                        ds = []
                        allPoints = BZ.fpoints
                        for vec in allPoints:
                            if areEqual(abs(dot(evec,vec)),norm(vec)): #axis and vec are parallel...don't want this one.
                                 ds.append(100)
                            else:
                                ds.append(norm(vec - evec*dot(evec,vec)))
                        allPoints = [point for (d,point) in sorted(zip(ds,allPoints),key = lambda x: x[0])]#sort by distance
                        pnto = allPoints[0]
                        #the plane to cut is in the plane of O and axis, but perpendiular to vector O.                   )
                        tempvec = cross(evec,pnto)
                        u1 = self.choose111(tempvec/norm(tempvec),eps)
                        BZ = self.cutCell(u1,0.0,BZ,eps)
                        pntp = dot(op,pnto)
                        tempvec = cross(evec,pntp)/norm(cross(evec,pntp))
                        if not allclose(tempvec,-u1):
                            u2 = self.choose111(tempvec/norm(tempvec),eps)
                            #print'\n\t2nd half of rot'
                            BZ = self.cutCell(u2,0.0,BZ,eps)
                    else: # -1: reflection/improper rotation
                        if len(where(areEqual(evals,-1.0))) > 1: evals = -evals #improper rotation
                        evec = evecs[:,where(areEqual(evals,-1.0))[0][0]]
                        u1 = self.choose111(evec,eps) 
                        BZ = self.cutCell(u1,0.0,BZ,eps)
#                     print 'iop',iop
#                     self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n'    
            elif areEqual(det(op),-1.0):
                inversion = True
        if inversion: #apply last of all
            if makesDups(array([[-1.,  0.,  0.], [ 0., -1.,  0.], [ 0.,  0., -1.]]),BZ):
                #can cut along any plane
                BZ = self.cutCell(array([1.0,0.0,0.0]),0.0,BZ,eps)
        self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n'
        hull = convexH(BZ.fpoints)
        BZ.volume = hull.volume
        self.IBZvolCut = det(self.B)/BZ.volume
        getBoundsFacets(BZ)
        BZ.fpoints = flatVecsList(BZ.facets) 
        BZ.center = sum(BZ.fpoints)/len(BZ.fpoints)
        print 'Vol BZ / Vol IBZ', self.IBZvolCut
        return BZ        

    def isAllInside(self,ds,eps):
        '''AllInside means on opposite side of the plane vs its normal vector, 
        and at least a distance of rpacking away from it'''
        allInside = zeros(len(ds),dtype = bool)
        for i,d in enumerate(ds):
            if self.method == 0:
                scale = sqrt(2.0)
            elif self.method == 1: 
                scale = sqrt(1.0)
            if d < - eps and abs(d)> (2*scale)*self.rpacking - eps: #point volume is all inside the true cell boundary
                allInside[i] = True
#         print 'Inside list', allInside
        return all(allInside)

    def isInsideExpanded(self,ds,eps):
        inside = zeros(len(ds),dtype = bool)
        for i,d in enumerate(ds):
            if d < - eps: #point is inside the expanded cell boundary
                inside[i] = True
#         print 'Exp in', inside
        return all(inside)

    def checkVertices(self,BZ,tMP,ds):
        nearPlanes = []
        nearDs = []
        for id, d in enumerate(ds):
            if abs(d)<=sqrt(2)*self.rpacking:
                for point in tMP.fpoints:
                    if dot(point,BZ.bounds[0][id]) < BZ.bounds[1][id]:
#                         print 'dot',dot(point,BZ.bounds[0][id]);sys.stdout.flush()
                        nearPlanes.append(id)
                        nearDs.append(d)
                        break
        return nearPlanes,nearDs
                        
    def facetsMathPrint(self,cell,label,axes = False, color = 'Red'):
        ''' Mathematica'''
        print '{}'.format(label),' = Graphics3D[{'+'{}'.format(color)+', Thick,{',
        for ifac, facet in enumerate(cell.facets):
            facet = list(trimSmall(array(facet)))
            print 'Line[{',
            for point in facet:
                print '{'+'{},{},{}'.format(point[0],point[1],point[2])+'},',
            print '{'+'{},{},{}'.format(facet[0][0],facet[0][1],facet[0][2])+'}', #back to first
            print '}]',
            if ifac < len(cell.facets)-1:
                print ',',
        print '}}',
        if axes:
            print ',Axes -> True,AxesLabel -> {"x", "y", "z"}]',
        else:
            print ']',   
        return    
            
    def facetsMeshMathPrint(self,cell):
        self.facetsMathPrint(cell,'p',True,'Red'); 
        print ';'
        print 'q=Graphics3D[{',
        cell.mesh = list(trimSmall(array(cell.mesh)))
        for ipoint,point in enumerate(cell.mesh):
            print 'Sphere[{' + '{},{},{}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(self.rpacking),
            if ipoint < len(cell.mesh) -1:
                print ',',
        print '}];',
    
    def facetsMeshVCMathPrint(self,BZ,MP):
        '''Put a mesh point VC at each mesh point'''
        self.facetsMathPrint(BZ,'s','True','Red'); #draw supecell voronoi cell
        showCommand = ';Show[s,'  
        print ';',
        for ipoint,point in enumerate(BZ.mesh):
            tCell = cell()
            tCell.facets = [[]]*len(MP.facets)
            for ifac,facet in enumerate(MP.facets):
                facet = list(trimSmall(array(facet)))
                temp = []
                for facetPoint in facet:
                    temp.append(facetPoint + point)
                tCell.facets[ifac] = temp
#             self.facetsMathPrint(tCell,'v{}'.format(ipoint),False,\
#                                  'discreteColors[{}][[{}]]'.format(len(BZ.mesh),ipoint));print ';',
            self.facetsMathPrint(tCell,'v{}'.format(ipoint),False,'RandomColor[]');print ';',

            showCommand += 'v{}'.format(ipoint)
            if ipoint < len(BZ.mesh)-1:
                showCommand += ','
#                 print ',', 
        showCommand += ']'
        print ';',
        print showCommand 
    
    def mathPrintPoints(self,points):
        '''As spheres'''
        print 'Graphics3D[{',
        for ipoint,point in enumerate(points):
            print 'Sphere[{' + '{},{},{}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(self.rpacking),
            if ipoint < len(points ) -1:
                print ',',
        print '}, Axes -> True, AxesLabel -> {"x", "y", "z"}]\n'
