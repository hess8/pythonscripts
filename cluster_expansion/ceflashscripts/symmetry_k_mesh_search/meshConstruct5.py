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
    orthdef, trimSmall, cosvecs,
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double,
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection,intoVoronoi,intoCell,
    reverseStructured,isInVoronoi,areParallel, among, addVec, getSGpointGroup)

def areEqual(x,y,eps):
    return abs(x-y)<eps

def addVec(vec,list,eps):
    '''adds a vector to a list of vectors if it's not in the list '''
    if among(vec,list,eps):
        return list
    else:
        list.append(vec)
    return list

def among(vec,list,eps):
    for lvec in list:
        if allclose(vec,lvec,atol=eps):
            return True
    return False 

def mathPrintPoints(points):
    print 's=Graphics3D[{',
    for ipoint,point in enumerate(points):
        print 'Sphere[{' + '{},{},{}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(0.05),
        if ipoint < len(points ) -1:
            print ',',
    print '}, Axes -> True, AxesLabel -> {"x", "y", "z"}]\n'

def mathPrintPlanes(planes):
    print 'r=RegionPlot3D[',
    for iplane, pvec in enumerate(planes):
        print '{}x+{}y+{}z<={}'.format(pvec[0],pvec[1],pvec[2],dot(pvec,pvec)), #write plane equations for mathematica
        if iplane < len(planes)-1:
            print '&&',
    print ',{x, -2, 2}, {y, -2, 2}, {z, -2, 2},PlotPoints -> 100]'

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
        if len(facet) == 3:
            return facet
        uvec = plane3pts(facet,eps)[0] #normal of facet, not cut plane.  These are the same only for bordersfacet
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
        return [point for (angle,point) in sorted(zip(angles,facet),key = lambda x: x[0])] 

def flatVecsList(vecsList,eps):
    '''Returns a flat list of vectors, from a list of lists that might have duplicates, 
    as in the case of facets on a polyhedron'''
    allpoints = []
    for isub,sub in enumerate(vecsList):
        for vec in sub:
            allpoints = addVec(vec,allpoints,eps)
    return allpoints

def plane3pts(points,eps):
    '''From the first 3 (noncollinear) points of a list of points, 
    returns a normal and closest distance to origin of the plane
    The closest distance is d = dot(r0,u).  The normal's direction 
    (degenerate vs multiplication by -1)
    is chosen to be that which points away from the side the origin is on.  If the 
    plane goes through the origin, then all the r's are in the plane, and no choice is made '''
    rcenter = sum(points)/len(points)
    r0 = points[0]; r1 = points[1]; r2 = points[2]
    vec = cross(r1-r0,r2-r0)
    nv = norm(vec)
    if nv < 5 * eps and len(points)>3: #two points are very close..use a 4th
        if norm(r1-r0)>norm(r2-r0): 
            r2 = points[3]
        else:
            r1 = points[3]
        vec = cross(r1-r0,r2-r0)
        nv = norm(vec)      
    u = vec/nv
    if dot(u,r0) < 0-eps:
        u = -u
    return u,dot(u,r0)

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

def newBounds(braggVecs,bndsLabels,grp,cell,eps):
    '''Find intersections planes:  newBraggs+current taken 2 at a time with e
    each newBragg, excluding duplicates in the triplet. 
    If any of these intersections are inside the current boundaries, then they 
    become new vertices, and the planes are added to boundaries. 
    
    If the new volume is the same as det(B), then we are done
    '''
    keepLabels = []
    checkNext = True
    nCurr = len(cell.bounds[0])
    allPlanes = [cell.bounds[0][i]*cell.bounds[1][i] for i in range(len(cell.bounds[0]))]
    allVerts = deepcopy(cell.fpoints)
    for ig  in grp:
        bndsLabels.append(ig)
    pairs = list(combinations(bndsLabels,2))
    iInt = 0
    for ig in grp:
        if ig ==5:
            'pause'
        for pair in pairs:
            planes3 = [braggVecs[ig]['vec']]
            if not (ig in pair):# or ig in keepLabels):
                planes3.append(braggVecs[pair[0]]['vec'])
                planes3.append(braggVecs[pair[1]]['vec'])
                intersPt = threePlaneIntersect(planes3)    
                if not intersPt is None:
                    iInt+=1
                    if not isOutside(intersPt,cell.bounds,eps)\
                        and not among(intersPt,allVerts,eps):
                            addVec(intersPt,allVerts,eps)
                            addVec(planes3[0],allPlanes,eps)
                            keepLabels.append(ig)
                            if pair[0]>=nCurr: 
                                keepLabels.append(pair[0])
                                addVec(planes3[1],allPlanes,eps)
                            if pair[1]>=nCurr: 
                                keepLabels.append(pair[1])
                                addVec(planes3[2],allPlanes,eps)           
    if len(allVerts)>0:
        #keep only the vertices that can be reached without crossing any plane
        newVerts = []
        for vert in allVerts:
            for plane in allPlanes:
                if dot(vert,plane)>dot(plane,plane)+eps:
                    break
            else: 
                newVerts.append(vert)
        #remove planes that don't host a vertex.
        newPlanes = []
        tryPlanes = deepcopy(allPlanes)
        for ip,plane in enumerate(tryPlanes):
            for vert in newVerts:
                if areEqual(dot(vert,plane),dot(plane,plane),eps):
                    addVec(plane,newPlanes,eps)
                    break
        cell.bounds = [[],[]]
        for plane in newPlanes:
            normPlane = norm(plane)
            cell.bounds[0].append(plane/normPlane)
            cell.bounds[1].append(normPlane)
        cell.fpoints = newVerts
#     mathPrintPlanes(allPlanes)
    if mod(len(cell.bounds[0]),2)!=0:
        sys.exit('Stop. Error in newBounds. BZ must have even number of planes')
    if len(cell.fpoints) >= 3:
#         mathPrintPoints(newVerts)
#         print 'Show[r,s]'
        print 'new Vol vs',convexH(newVerts).volume,cell.volume
        checkNext = not areEqual(convexH(newVerts).volume,cell.volume,eps/4.0) #should be some power greater than 1
    else:
        checkNext = True
    return checkNext,bndsLabels,cell

def getInterscPoints(planes,eps):
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
            interscPoints = addVec(interscP,interscPoints,eps)
    for ipoint, vec in enumerate(interscPoints):
        interscPmags.append(norm(vec))
    interscPoints = [point for (mag,point) in sorted(zip(interscPmags,interscPoints),key = lambda x: x[0])] #sorted by distance from origin 
    return interscPoints

def getFacetsPoints(cell,eps):
    cell.facets = [[]]*len(cell.bounds[0])
#     facetsPoints = []
#     for point in interscPoints:
#         if not isOutside(point,cell.bounds,eps):
#             facetsPoints.append(point)
    #arrange intersection points in each plane according to their angular order in the plane
    for iplane, uvec in enumerate(cell.bounds[0]):
        pvec = uvec*cell.bounds[1][iplane]
        facetvecs = []
        for i, vec in  enumerate(cell.fpoints):
            if onPlane(vec,pvec,eps) and not isOutside(vec,cell.bounds,eps):
                facetvecs = addVec(vec,facetvecs,eps)          
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
                    
def makesDups(op,cell,eps):
    '''Applies symmetry operator to all facet points. If any facet points are 
    moved on top of other facet points, then return true'''
    points = cell.fpoints
    for i in range(len(points)):
        rpoint = dot(op,points[i])
        if allclose(points[i],rpoint,atol=eps):
            break #no rotation effect
        otherLabels = range(len(points))
        otherLabels.pop(i)
        for j in otherLabels:
            if allclose(rpoint,points[j],atol=eps):
                return True
    return False

def getVorCell(LVs,cell,eps):
    '''Boundaries and vertices of Voronoi cell'''
    braggVecs = getBraggVecs(LVs)
    
    igroup = 1
#     mathPrintPoints(braggVecs[:]['vec'])
    checkNext = True
    gstart,ng = magGroup(braggVecs,1,eps) # group of smallest bragg plane vectors
    boundsLabels = range(ng)
    for i in boundsLabels:
        vec = braggVecs[i]['vec']; mag = norm(vec)
        cell.bounds[0].append(vec/mag); cell.bounds[1].append(mag)
    #get any intersections between these planes
    cell.fpoints = getInterscPoints(cell.bounds,eps)
    while checkNext:
        igroup += 1
#         print 'igroup',igroup
        gstart,ng = magGroup(braggVecs,igroup,eps)
        nextGroup = range(gstart,gstart+ng)
#         checkNext,cell = newBoundsifInside(braggVecs,nextGroup,cell,eps)
        checkNext,boundsLabels,cell = newBounds(braggVecs,boundsLabels,nextGroup,cell,eps)
        if igroup > 10:
            sys.exit('Stop. getVorCell igroup >10.  This seems to be an infinite while loop.  Try reducing eps ')
#     interscPoints = getInterscPoints(cell.bounds,eps)
    #check to see if an
    #write plane equations for Mathematica:
#         for iplane, uvec in enumerate(cell.bounds[0]):
#             pvec = uvec*cell.bounds[1][iplane]
#             print '{}x+{}y+{}z<={}&&'.format(pvec[0],pvec[1],pvec[2],dot(pvec,pvec)) #write plane equations for mathematica
    cell = getFacetsPoints(cell,eps)
    cell.fpoints = flatVecsList(cell.facets,eps)
    cell.center = sum(cell.fpoints)/len(cell.fpoints)
    cell.volume = convexH(cell.fpoints).volume
#     for point in cell.fpoints:
#         print 'fpoint {:20.17f} {:20.17f} {:20.17f}'.format(point[0],point[1],point[2])  
    return cell
             
def getBraggVecs(LVs):
    '''The Bragg vectors are halfway from the origin to a lattice point.
    The planes normal to some of these will be bounding planes of the Voronoi cell '''
    braggVecs = zeros(5*5*5-1,dtype = [('vec', '3float'),('mag', 'float'),('dep', 'S15')])  
    #Exclude (0,0,0) in array dimensions (-1)
    ipoint = 0
    for i in range(-2,3):
        for j in range(-2,3):
            for k in range(-2,3):
                if not (i==j==k==0):
                    vec = trimSmall(0.5*(i*LVs[:,0] + j*LVs[:,1] + k*LVs[:,2]))
#                     vec = trimSmall(0.5*dot(LVs,array([i,j,k])))
                    
                    braggVecs[ipoint]['vec'] = vec
                    braggVecs[ipoint]['dep'] = '{},{},{}'.format(i,j,k)
                    braggVecs[ipoint]['mag'] = norm(vec)
                    ipoint+=1
    braggVecs.sort(order = 'mag')
#     for vec in braggVecs['vec']:
#         print 'bragg vector', vec
    return braggVecs

# def intsPlLinSeg(u,ro,r1,r2,eps):
#     '''Intersection between a plane through the origin and a line.
#     A plane through the origin is given by dot(r,u) = ro.  
#     A line segment between r1 and r2 is given by vectors r = r1+t(r2-r1), for t in (0,1) 
#     Combining these:  r1u = dot(r1,u).  r2u = dot(r2,u).  Then t = (ro-r1u)/(r2u-r1u).
#     So if t is in (0,1), then we have an intersection'''
#     r1u = dot(r1,u)
#     r2u = dot(r2,u)
#     if areEqual(r1u,r2u,eps):
#         return False, None
#     else:
#         t = (ro-r1u)/(r2u-r1u)
#         if 0 + eps < t < 1.0-eps:
#             return True, trimSmall(r1 + t*(r2-r1)) 
#         else:
#             return False, None    

def intsPlLinSeg(u,ro,r1,r2,eps):
    '''Intersection between a plane and a line.
    A plane is given by dot(r,u) = ro.  
    A line segment between r1 and r2 is given by vectors r = r1+t(r2-r1), for t in (0,1) 
    Combining these:  r1u = dot(r1,u).  r2u = dot(r2,u).  Then t = (ro-r1u)/(r2u-r1u).
    So if t is in (0,1), then we have an intersection'''
    r1u = dot(r1,u)
    r2u = dot(r2,u)
#     if areEqual(r1u,r2u,eps):
#         return False, None
    t = (ro-r1u)/(r2u-r1u)
    return r1 + t*(r2-r1) 

def getBoundsFacets(cell,eps,rpacking = None): 
    '''After cutting, we have new facets, and the facet planes are the boundaries'''
    cell.bounds = [[],[]] 
    
    for ifac, facet in enumerate(cell.facets):
        facCenter = sum(facet)/len(facet)
        u,ro = plane3pts(facet,eps)
        if areEqual(ro,0,eps) and dot(facCenter - cell.center, u )< 0-eps: #plane goes through origin...choose a normal that points "outside" the BZ
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
#     def __init__(self):
#         self.BZboundaries = [[],[]] #for BZ Voronoi cell or IBZ.  written as [u's] , [ro's]
#         self.MPboundaries = [[],[]] #for mesh point Voronoi cell
#         cell.facets = [] #facet points for BZ Voronoi cell or IBZ
#         self.MPfacets = [] #facet points for mesh point Voronoi cell
        
    def meshSym(self,A,B,totatoms,aTypes,postype,aPos,targetNmesh,meshtype,path,method):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        self.B = B
        print '\nB (Recip lattice vectors as columns',B
        print 'method',method
        [self.symops,self.nops] = getGroup(self.B)
#         for iop in range(self.nops):
#             print 'iop',iop
#             for i in range(3):
#                 for j in range(3):
#                     print '{:20.17f}'.format(self.symops[i,j,iop])
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
        self.path = path
        self.method = method
            #0: exact: use vertices of mesh voronoi cell that are closest/farthest 
            #         from the IBZ center origin to check if the point's volume is cut. 
            #         Cut the VC to determine the volume contribution      
            #0.5 approx  If cut volume is less than 50%, distribute weight to neighbors of equivalent points
            #         If point is outside of first BZ, then translate by G vector to get it inside.  Points inside but not in IBZ
            #        use point symmetries to get in.  This applies to kpoints in IBZ but near corners as well. 
        print 'Number of desired points:', targetNmesh
        print 'Symmetry operations:', self.nops
        vol = abs(det(B))
        self.ravg = (vol/targetNmesh)**(1/3.0)
        self.edgeFactor = 3.0
        self.rmin = 0.8*self.ravg #
        self.rmax = self.edgeFactor*self.ravg #cutoff for forces, neighbors in other cells. 
        eps = self.ravg/2000
        BZ = cell() #instance
        BZ.volume = vol
        BZ = getVorCell(self.B,BZ,eps)
#         print 'Vornonoi cell'; self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n'
        self.vorCell = BZ
#         self.facetsMathPrint(BZ,'p',True,'Red') 
#         IBZ = self.getIBZ(BZ,eps) #now irreducible BZ  
        print 'Temporary: no symmetry reduction'
        IBZ = self.vorCell #temp code!!! no symmetry reduction
        self.IBZvolCut = 1.0 #temp code!!! no symmetry reduction
        self.meshCubic(IBZ,meshtype,eps) #cub, fcc, bcc   
#         self.meshCubic(IBZ,'fcc',eps)
#         self.meshCubic(IBZ,'cub',eps)
#         self.triFaces()
#         self.meshCubic('bcc')

#         self.meshCubic('cub')   
        self.writeKpoints(IBZ)       
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
                         .format(kpoint[0],kpoint[1],kpoint[2],cell.weights[ik]*1e6))
        self.writefile(lines,'KPOINTS')         
   
    def meshCubic(self,IBZ,type,eps):
        '''Add a cubic mesh to the interior, . If any 2 or 3 of the facet planes are 
        orthogonal, align the cubic mesh with their normals.       
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
'''
        
        a = 1.0
        cubicLVs = identity(3)
        #test facet points for orthogonality
        rs = []
        pairs = []
        triples = []
        uvecs = IBZ.bounds[0]
        for i in range(len(uvecs)):
            if areEqual(norm(uvecs[i]),0.0,eps): break
            rs.append(norm(uvecs[i]))
            for j in range(i,len(uvecs)):
                if areEqual(norm(uvecs[j]),0.0,eps): break
                if areEqual(dot(uvecs[i],uvecs[j]),0.0,eps):
                    pairs.append([uvecs[i],uvecs[j]])

        for ip,pair in enumerate(pairs):
            for i in range(len(uvecs)):
                if areEqual(norm(uvecs[i]),0.0,eps): break
                if areEqual(dot(pair[0],uvecs[i]),0.0,eps) and areEqual(dot(pair[1],uvecs[i]),0.0,eps):
                    triple = deepcopy(pair)
                    triple.append(uvecs[i])
                    triples.append(triple)
                    break
        #Define basis vectors for cubic lattice:
        Lsum= [] #length of vectors in pair or triplet
        if len(triples)>0:
            print 'At least one triplet of orthogonal plane normals found:',triples[0]
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
            print 'At least one pair of orthogonal plane normals found:', pairs[0]
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
            print 'no orthogonal plane normals pairs found.'
        if type == 'fcc':
            volKcubConv = det(self.B)/self.nTarget*4
            aKcubConv = volKcubConv**(1/3.0)
            cubicLVs = cubicLVs * aKcubConv
            sites = [array([0, 0 , 0]), 1/2.0*(cubicLVs[:,1]+cubicLVs[:,2]),\
                     1/2.0*(cubicLVs[:,0]+cubicLVs[:,2]), 1/2.0*(cubicLVs[:,0]+cubicLVs[:,1])]
            primLVs = transpose(array(sites[1:]))
            self.rpacking = 1/2.0/sqrt(2)*aKcubConv
            pf = 4*4/3.0*pi*(1/2.0/sqrt(2))**3  #0.74
#             fccFacets = 
        elif type == 'bcc':
            volKcubConv = det(self.B)/self.nTarget*2
            aKcubConv = volKcubConv**(1/3.0)
            cubicLVs = cubicLVs * aKcubConv
            sites = [array([0, 0 , 0]), 1/2.0*(cubicLVs[:,0]+cubicLVs[:,1]+cubicLVs[:,2])]
            primLVs = transpose(array([array(-1/2.0*(cubicLVs[:,0]+cubicLVs[:,1]+cubicLVs[:,2])),\
                        array(1/2.0*(cubicLVs[:,0]-cubicLVs[:,1]+cubicLVs[:,2])),\
                        array(1/2.0*(cubicLVs[:,0]+cubicLVs[:,1]-cubicLVs[:,2]))]))
            self.rpacking = sqrt(3)/4.0*aKcubConv
            pf = 2*4/3.0*pi*(sqrt(3)/4.0)**3 #0.68
        elif type == 'cub':
            volKcubConv = det(self.B)/self.nTarget
            aKcubConv = volKcubConv**(1/3.0)
            cubicLVs = cubicLVs * aKcubConv
            sites = [array([0, 0 , 0])]
            primLVs = cubicLVs
            self.rpacking = aKcubConv/2
            pf = 4/3.0*pi*(1/2.0)**3 #0.52
        else:
            sys.exit('stop. Type error in meshCubic.')
#         primLVs = primLVs * aKcubConv
#         cubicLVs = cubicLVs * aKcubConv
#         sites = [site * aKcubConv for site in sites]
#         self.rpacking = self.rpacking*aKcubConv
        self.Vsphere = 4/3.0*pi*self.rpacking**3
        print 'rpacking',self.rpacking
        IBZ = getBoundsFacets(IBZ,eps,self.rpacking) #adds expanded cell
        MP = cell()
        MP.volume = volKcubConv/len(sites)
        MP = getVorCell(primLVs,MP,eps)
        MP.volume = convexH(MP.fpoints).volume
#         print 'Mesh Point Voronoi cell',MP.facets; self.facetsMathPrint(MP,'p',True,'Red');print ';Show[p]\n'        
        #Find the extremes in each cubLV direction:
        intMaxs = [] #factors of aKcubConv
        intMins = []
        for i in range(3):
#             print 'cubic',i,cubicLVs[:,i]
            projs = []
            for point in IBZ.fpoints:
                projs.append(dot(cubicLVs[:,i],point)/aKcubConv**2)
            intMaxs.append(int(ceil(max(projs)))+1)
            intMins.append(int(floor(min(projs))))
        print 'Maxes',intMaxs
        print 'Mins',intMins       
        #Create the cubic mesh inside the irreducible BZ
        IBZ.mesh = []
        IBZ.weights = []
        kptsRed = []
        wgtsRed = []
        dsRed = []
        nRed = 0
        weightsInside = 0
        nInside = 0
        weightsOnePlane = 0
        nOnePlane = 0
        weightsCuts = 0
        nCut = 0
        nDummy = 0
#         offset = array([0.5,0.5,0.5])*aKcubConv
#         print 'bounds'
#         for i in range(len(IBZ.bounds[0])):
#             print IBZ.bounds[0][i],IBZ.bounds[1][i]               
        ik = 0 
        #begin MP facets printing
#         self.facetsMathPrint(IBZ,'s','True','Red'); print ';', #draw supecell voronoi cell before loop
        showCommand = 'Show[s,' 
        #done with start of MP facets printing
        
#         shift = aKcubConv*array([0,0,0])/2.0
        shift = aKcubConv*array([1,1,1])/2.0
        for i in range(intMins[0],intMaxs[0]):
            for j in range(intMins[1],intMaxs[1]):
                for k in range(intMins[2],intMaxs[2]):
                    lvec = i*cubicLVs[:,0]+j*cubicLVs[:,1]+k*cubicLVs[:,2]
#                     if i==j==k==-1:
#                         'pause'
                    for iS, site in enumerate(sites):
                        ik+=1
                        kpoint = lvec + site + shift
#                         print 'ik',ik
#                         print 'test',[i,j,k],iS,kpoint
#                         ds = self.dToPlanes(kpoint,IBZ.expBounds)

#                         print 'kpoint',ik,i,k,j,kpoint

#                         print 'ds',ds;print
#                         inExpanded,centerInside,allInside = self.boundStatus(ds,eps)
                        cutMP = self.prepCutMP(MP,kpoint)  #cut MP is displaced by kpoint from MP
                        self.dotsVsRos = zeros(len(MP.fpoints),len(IBZ.bounds[0]))#Where each fpoint lies vs the plane: 0 on the plane >0 outside
                        inExpanded,centerInside,allInside = self.boundStatusExact(cutMP,IBZ,eps)
#                         print 'inExpanded,centerInside,allInside', inExpanded,centerInside,allInside
                        if allInside:
#                             print 'allInside'
                            IBZ.mesh.append(kpoint)
                            IBZ.weights.append(self.IBZvolCut)
                            weightsInside += self.IBZvolCut
                            nInside += 1
                        elif inExpanded:
#                             print 'inExpanded'
#                             if centerInside: print 'centerInside'
                            #change to d's from real cell, the borders of the IBZ
#                             dsIBZ = [d + sqrt(2)*self.rpacking for d in ds]
#                             print '\n\nkpoint',kpoint, ik, [i,j,k]; sys.stdout.flush()
#                             print 'ds',dsIBZ                           
#                             cutMP = self.prepCutMP(MP,kpoint)  #cut MP is displaced by kpoint from MP
#                             print;print 'ik',ik
#                             if ik == 115:
#                                 'pause'
                            for iplane, uvec in enumerate(IBZ.bounds[0]):
                                ro = IBZ.bounds[1][iplane]                                      
                                cutMP = self.cutCell(uvec,ro,cutMP,eps) # we always keep the part that is "inside", opposite u
#                                 self.facetsMathPrint(cutMP,'p',True,'Red'); print ';Show[p]\n' 
                                if len(cutMP.facets) <4:
                                    cutMP.volume = 0.0
                                if cutMP.volume == 0.0: #(outside IBZ. happens in oblique corners of expanded cell)
                                        break
                            if len(cutMP.facets)>=4:
                                cutMP.volume = convexH(cutMP.fpoints).volume
#                             print 'volume cut',cutMP.volume, cutMP.volume > eps**3
                            if cutMP.volume > eps**3:
                                weight = self.IBZvolCut*cutMP.volume/MP.volume
                                if self.method > 0  and not centerInside:  #kpoints outside of the IBZ have their weight redistributed.  Not necessary 0.0 < cutMP.volume < 0.5*MP.volume
                                    nDummy +=1
                                    nRed += 1
#                                     kptsRed.append(kpoint)
                                    kptsRed.append(cutMP.center) #If the mech VC is cut, may as well use its center.  
                                    wgtsRed.append(weight)
                                    dsRed.append(ds)
#                                     print 'kpoint',ik,iS,'cut',nDummy,i,k,j,kpoint,'vol',weight*MP.volume
                                else:   
#                                     IBZ.mesh.append(kpoint)
                                    IBZ.mesh.append(cutMP.center) #If the mesh VC is cut, may as well use its center.  
                                    IBZ.weights.append(weight)   
                                    weightsCuts += weight
                                    nCut += 1
                                    nDummy +=1
#                                     print 'kpoint',ik,iS,'ndum',nDummy,i,k,j,kpoint,'vol',weight*MP.volume
                                #####MP facet printing loop line                                
#                                 showCommand = self.cutMPCellMathPrint(IBZ,cutMP,kpoint,ik,showCommand)
                                #####end MP facet printing loop entry
        if showCommand[-1] == ',': showCommand = showCommand[:-1]
        showCommand += ']' 
        print ';', 
        print showCommand 
        #end MP facets printing
        #redistribute or reassign low-volume weights
        if  len(kptsRed) <= len(IBZ.mesh) and len(kptsRed)>0:
            weightRed = sum(wgtsRed)
            print 'Weights redistributed', nRed,weightRed, 'Average per point', weightRed/float(nRed)        
            IBZ = self.redistribWgt(kptsRed,wgtsRed,dsRed,IBZ,eps)
        elif len(kptsRed)>0:
            print 'Too few inside kpoints to redistribute: keeping all border points'
            IBZ = self.addToKpts(kptsRed,wgtsRed,IBZ)
            nCut += len(kptsRed)
            weightsCuts += sum(wgtsRed)
        print 'Weights allinside pts', nInside, weightsInside
#         print 'Volume inside pts', weightsInside*MP.volume
        if nOnePlane>0: print 'Weights one plane',nOnePlane,weightsOnePlane, weightsOnePlane/float(nOnePlane)
        print 'Weights in IBZ cuts', nCut,weightsCuts, 'Average per cut VC', weightsCuts/float(nCut)
        print 'Total volume in weights:',  sum(IBZ.weights)*MP.volume, 'from ', (nCut + nOnePlane + nInside),'points'
        print 'BZ volume:', det(self.B),'\n'
        self.facetsMeshMathPrint(IBZ); print ';Show[p,q]\n'
        self.facetsMeshVCMathPrint(IBZ,MP)
        return
    
    def redistribWgt(self,kptsRed,wgtsRed,dsRed,IBZ,eps):
        ''' 0.5 approx  If the kpoint is outside, distribute the portion of its Vornoi cell weight that is inside to neighbors of equivalent points
            if point is outside of Voronoi cell (first BZ), then translate by G vector to get it inside.  Points inside the Voronoi cell but not in IBZ
            use point symmetries to get in the IBZ.   
         '''  
        for ik,kpoint in enumerate(kptsRed):
#             print;print 'kpoint',ik, kpoint
            dists = []
            neighWgts = []
            ds = dsRed[ik]
            if ik == 4:
                'pause'
            neighs, neighLbls = self.getNeighbors(kpoint,ds,IBZ,eps)         
#             print 'neighs',neighs,neighLbls
            totNsWeights = 0
            onTop = False 
            for iN, neigh in enumerate(neighs):
                d = norm(kpoint-neigh)
                if areEqual(d,0.0,eps): #kpoint lies on top of a neighbor which gets all the weight
                    IBZ.weights[neighLbls[iN]] += wgtsRed[ik]
                    onTop = True
                    break
                dists.append(d)
                neighWgt = 1/d
                neighWgts.append(neighWgt)
                totNsWeights += neighWgt
            if not onTop:
                for iN, neigh in enumerate(neighs): 
                    IBZ.weights[neighLbls[iN]] += wgtsRed[ik] * neighWgts[iN]/totNsWeights
        return IBZ
    
    def intoIBZ(self,kpoint,ds,IBZ,eps):   
        '''For a point inside the Voronoi cell, use point symmetry to get into IBZ.
        If a point is just outside the Voronoi cell, then translate it by -2*(plane vector) of the bragg plane
        that it is closest to. Then use point symmetry to get into IBZ. This should turn out to be
        equivalent to a reflection about the plane, or a reflection and translation along the plane.
        
        At least one bragg plane must be part of the IBZ bounds.  These do not pass through
        the origin'''
        
        sympoints = []
        if isOutside(kpoint,self.vorCell.bounds,eps):
            kpoint = intoVoronoi(kpoint,self.B)
        for iop in range(self.nops):
            op = self.symops[:,:,iop]
            kpoint2 = dot(op,kpoint)
            sympoints.append(kpoint2)
            if not isOutside(kpoint2,IBZ.bounds,eps):
                return kpoint2
        else:
            self.mathPrintPoints(sympoints)
            sys.exit("Stop. intoIBZ: symm ops don't return a kpoint inside the IBZ")
                    
    def getNeighbors(self,kpoint,ds,IBZ,eps):
        '''This is for a kpoint just outside the IBZ. 
        Search a sphere around the kpoint and collect neighbors.
        Then if outside the IBZ, move point into IBZ by symmetry.  Search another sphere.
        Return the closest three neighbors
        '''
        neighR = 3.0*self.rpacking
        kpSymm = self.intoIBZ(kpoint,ds,IBZ,eps)
        neighs,neighLbls = self.searchSpheres([kpoint,kpSymm],IBZ,neighR)
        return neighs,neighLbls 
    
    def searchSpheres(self,klist,IBZ,R):
        neighs = []
        neighLbls = []
        for ip,meshPt in enumerate(IBZ.mesh):
            for kpoint in klist:
                if norm(meshPt-kpoint) <= R:
                    neighs.append(meshPt)
                    neighLbls.append(ip)
                    break
        return neighs,neighLbls
        
    def addToKpts(self,kptsRed,wgtsRed,IBZ):
        for ik,kpoint in enumerate(kptsRed):
            IBZ.mesh.append(kpoint)
            IBZ.weights.append(wgtsRed[ik])
        return IBZ
        
    def prepCutMP(self,MP,kpoint):
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
        ftemp = deepcopy(cell.facets)       
#         print;print 'u',u,ro
#         if u[0] <0 and u[1]<0 and u[2]>0:
#             print;print 'u',u,ro
#         if allclose(u,array([-0.57735042,  0.21442248 ,-0.7878385]),atol=eps):
#             'pause'
        for ifac, facet in enumerate(cell.facets):
#             print 'ifac',ifac
#             if ifac == 7:
#                 'pause'
            newfacet = []
            #first, determine if the plane intersects any facet points or 
            #line segments, by checking the change of sign of dot products
            signs = []
            for point in facet:
                pu = dot(u,point)
                if pu > ro + eps:
                    signs.append(1.0)
                elif areEqual(pu,ro,eps):
                    signs.append(0.0)
                else:
                    signs.append(-1.0)
            if -1 in signs and 1 in signs: #facet is cut
                for ip,pointi in enumerate(facet):
                    if signs[ip] == 0.:# then this point is on the cut plane 
                        newfacet.append(pointi)
                        bordersFacet = addVec(pointi,bordersFacet,eps)
                    elif signs[ip] == -1.0:
                        newfacet.append(pointi)
                    else: #outside
                        allRemoved = addVec(pointi,allRemoved,eps)
                    jp = ip + 1
                    if jp == len(facet): jp = 0 
                    if signs[ip]!=0 and signs[ip] == -signs[jp]: #will have intersection between point1 and next point on facet
                        pointj = facet[jp]
#                         [intersect, rinters] = intsPlLinSeg(u,ro,pointi,pointj,eps) #         
#                         if intersect:
#                             newfacet.append(rinters)
#                             bordersFacet = addVec(rinters,bordersFacet,eps)
#                         else:
#                             print'error in cutCell'
                        interscP = intsPlLinSeg(u,ro,pointi,pointj,eps) #         
                        newfacet.append(interscP)
                        bordersFacet = addVec(interscP,bordersFacet,eps)

                if len(newfacet) >= 3:
#                     print 'newfacet',newfacet
                    ftemp[ifac] = orderAngle(newfacet,eps)
            else: #mark for removal all points that are outside of the plane
                for i, sgn in enumerate(signs):
                    if sgn == 1.0:
                        allRemoved = addVec(facet[i],allRemoved,eps)
                    elif sgn == 0.0:
                        bordersFacet = addVec(facet[i],bordersFacet,eps)         
        if len(allRemoved)>0:
            ftemp2 = deepcopy(ftemp)
            for i2, facet in enumerate(ftemp):
                nextbreak = False
                for ip,point in enumerate(facet):
                    for rpoint in allRemoved: 
                        if allclose(point,rpoint,atol=eps):
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
            points = flatVecsList(cell.facets,eps)
            for point in points:
                if areEqual(dot(u,point),ro,eps):
                    bordersFacet = addVec(point,bordersFacet,eps)            
            #Order by angle in the facet
            if len(bordersFacet)> 2:
#                 print 'bordersfacet'
                ftemp.append(orderAngle(bordersFacet,eps))
            cell.facets = ftemp
            cell.fpoints = flatVecsList(cell.facets,eps) 
            if len(cell.fpoints)== 0:
                cell.volume = 0
            else:
                cell.center = sum(cell.fpoints)/len(cell.fpoints)
#                 print'center',cell.center
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
        has complex eigenvalues and a determinate .... 
           
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
        
        if not areEqual(det(self.B),BZ.volume,eps):
            subprocess.call(['echo', '\nError: Voronoi cell volume {} is not equal to the parallelepiped volume {}'.format(BZ.volume,det(self.B))])
            sys.exit('Stop.')
        print '\n\nReducing Brillouin zone by symmetry'
#         self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n' 
         
        inversion = False
        for iop in range(self.nops):
            op = self.symops[:,:,iop] 
#             print 'symop',iop;print op ;print          
            if abs(trace(op))< 3.0-eps: #skip E and inverse
                evals,evecs = eig(op)
                evecs = array([evec for evec in evecs])
                if areEqual(det(op),-1.0,eps) and not allclose(imag(evecs),zeros((3,3)),atol=eps):
#                   Change improper rotation to proper one'
                    op = -op
                    evals,evecs = eig(op)
                    evecs = array([evec for evec in evecs])
#                 if areEqual(det(op),1.0)  : #rotation
#                     print 'Rotation',  real(evecs[:,where(areEqual(evals,1.0))[0][0]])
#                 else:
#                     print 'Reflection',real(evecs[:,where(areEqual(evals,-1.0))[0][0]])
                if makesDups(op,BZ,eps): #does this operation cause current facet points to move to other facet points
                    if areEqual(det(op),1.0,eps)  : #rotation
                        evec = evecs[:,where(areEqual(evals,1.0,eps))[0][0]] #axis
                        ipt = 0
                        #choose a facet point that is close to the rotation axis to avoid unusual cutting planes
                        ds = []
                        allPoints = BZ.fpoints
                        for vec in allPoints:
                            if areEqual(abs(dot(evec,vec)),norm(vec),eps): #axis and vec are parallel...don't want this one.
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
                        tempvec = cross(evec,pntp)/norm(cross(evec,pntp))#2nd cut plane for roation
                        if not allclose(tempvec,-u1,atol=eps): #don't cut again if this is a Pi rotation
                            u2 = self.choose111(tempvec/norm(tempvec),eps)
                            BZ = self.cutCell(u2,0.0,BZ,eps)
                    else: # -1: reflection/improper rotation
                        if len(where(areEqual(evals,-1.0,eps)) )> 1: evals = -evals #improper rotation
                        evec = evecs[:,where(areEqual(evals,-1.0,eps))[0][0]]
                        u1 = self.choose111(evec,eps) 
                        BZ = self.cutCell(u1,0.0,BZ,eps)
#                     self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n'   
            elif areEqual(det(op),-1.0,eps):
                inversion = True
        if inversion: #apply last of all
            if makesDups(array([[-1.,  0.,  0.], [ 0., -1.,  0.], [ 0.,  0., -1.]]),BZ,eps):
                #can cut along any plane
                BZ = self.cutCell(array([1.0,0.0,0.0]),0.0,BZ,eps)
#         self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n'
        BZ.volume = convexH(BZ.fpoints).volume
        self.IBZvolCut = det(self.B)/BZ.volume
        getBoundsFacets(BZ,eps)
        BZ.fpoints = flatVecsList(BZ.facets,eps)
        BZ.center = sum(BZ.fpoints)/len(BZ.fpoints)
        print 'Vol BZ / Vol IBZ', self.IBZvolCut
        return BZ

    def boundStatusExact(self,MP,IBZ,eps):
        '''Check all the facets vs the IBZ bounds'''        
#         inExpanded = zeros(nbounds,dtype = bool)
#         centerInside = zeros(nbounds,dtype = bool)
#         allInside = zeros(nbounds,dtype = bool)
        for ipoint,point in enumerate(MP.fpoints):
            for iu,uvec in enumerate(IBZ.bounds[0]):
                ro = IBZ.bounds[0][iu]
                pu = dot(u,point)
                self.dotsVsRos[ipoint,iu] = pu - ro
#                     if pu > ro + eps:
#                         dots.append(1.0)
#                     elif areEqual(pu,ro,eps):
#                         dots.append(0.0)
#                     else:
#                         dots.append(-1.0)
#                     maxsign = max(signs); minsign = min(signs)
#                     if maxsign == 1.0:
#                 if max(dots)
#                     allInside[i] = True 
#                     
#                 inExpanded[i] = True
#                 if abs(d) >= (2*scale)*self.rpacking - eps: #point volume is all inside the true cell boundary
#                     allInside[i] = True 
#                     centerInside[i] = True 
#                 elif abs(d) >= (1*scale)*self.rpacking - eps:
#                     centerInside[i] = True
        mindots = amin(dotsVsRos)
        maxdots = amax(dotsVsRos)
        toCut = False; centerInside = False; allInside = False
        if mindots >= 0-eps: #all of MP is outside of IBZ or at most touching boundary
            return 'outside'
        elif maxdots > 0 + eps and mindots < 0 -eps:
            return 'cut'
        elif maxdots < 0 - eps:
            return 'allInside' 

#     def boundStatus(self,ds,eps):
#         '''The d's are from the expanded boundary
#         allInside: the center is on the opposite side of the plane vs its normal vector, 
#         and at least a distance of rpacking*scale away from it.  
#         centerInside: the center is inside or on the boundaries
#         inExpanded: the center is inside the expanded boundaries'''
#         inExpanded = zeros(len(ds),dtype = bool)
#         centerInside = zeros(len(ds),dtype = bool)
#         allInside = zeros(len(ds),dtype = bool)
#         for i,d in enumerate(ds):
#             if self.method < 1:
#                 scale = sqrt(2.0) #finding if a voronoi cell edge is inside.  Does this work for bcc? 
#             elif self.method == 1: 
#                 scale = sqrt(1.0)
#             if d < -eps:
#                 print 'd',abs(d),(1*scale)*self.rpacking - eps,(2*scale)*self.rpacking - eps
#                 inExpanded[i] = True
#                 if abs(d) >= (2*scale)*self.rpacking - eps: #point volume is all inside the true cell boundary
#                     allInside[i] = True 
#                     centerInside[i] = True 
#                 elif abs(d) >= (1*scale)*self.rpacking - eps:
#                     centerInside[i] = True
#         return all(inExpanded), all(centerInside), all(allInside)        
#         

#     def isAllInside(self,ds,eps):
#         '''AllInside means on opposite side of the plane vs its normal vector, 
#         and at least a distance of rpacking away from it'''
#         allInside = zeros(len(ds),dtype = bool)
#         for i,d in enumerate(ds):
#             if self.method == 0:
#                 scale = sqrt(2.0)
#             elif self.method == 1: 
#                 scale = sqrt(1.0)
#             if d < - eps and abs(d)> (2*scale)*self.rpacking - eps: #point volume is all inside the true cell boundary
#                 allInside[i] = True
#         return all(allInside)

#     def isInsideExpanded(self,ds,eps):
#         inside = zeros(len(ds),dtype = bool)
#         for i,d in enumerate(ds):
#             if d < - eps: #point is inside the expanded cell boundary
#                 inside[i] = True
#         return all(inside)

    def checkVertices(self,BZ,tMP,ds):
        nearPlanes = []
        nearDs = []
        for id, d in enumerate(ds):
            if abs(d)<=sqrt(2)*self.rpacking:
                for point in tMP.fpoints:
                    if dot(point,BZ.bounds[0][id]) < BZ.bounds[1][id]:
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
    
    def facetsMathToStr(self,strOut,cell,label,axes = False, color = 'Red'):
        ''' Mathematica'''
        strOut += '{}'.format(label)+' = Graphics3D[{'+'{}'.format(color)+', Thick,{'
        for ifac, facet in enumerate(cell.facets):
            facet = list(trimSmall(array(facet)))
            strOut += 'Line[{'
            for point in facet:
                strOut += '{'+'{},{},{}'.format(point[0],point[1],point[2])+'},'
            strOut += '{'+'{},{},{}'.format(facet[0][0],facet[0][1],facet[0][2])+'}' #back to first
            strOut += '}]'
            if ifac < len(cell.facets)-1:
                strOut += ','
        strOut += '}}'
        if axes:
            strOut += ',Axes -> True,AxesLabel -> {"x", "y", "z"}]'
        else:
            strOut += ']'  
        return strOut
            
    def facetsMeshMathPrint(self,cell):
        '''Writes spheres for Mathematica drawing'''
        self.facetsMathPrint(cell,'p',True,'Red'); 
        print ';'
        print 'q=Graphics3D[{',
        cell.mesh = list(trimSmall(array(cell.mesh)))
        for ipoint,point in enumerate(cell.mesh):
            print 'Opacity[0.7],Sphere[{' + '{},{},{}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(self.rpacking),
            if ipoint < len(cell.mesh) -1:
                print ',',
        print '}];',
   
    def facetsMeshVCMathPrint(self,BZ,MP):
        '''Put a mesh point VC at each mesh point'''
        mathOut = open('mathOut','w') #both prints and writes string to file 
        strOut = ''
        self.facetsMathPrint(BZ,'s','True','Red'); #draw supecell voronoi cell
        strOut = self.facetsMathToStr(strOut,BZ,'s','True','Red'); #draw supecell voronoi cell
        showCommand = ';Show[s,'  
        print ';',;strOut+=';'
        for ipoint,point in enumerate(BZ.mesh):
            tCell = cell()
            tCell.facets = [[]]*len(MP.facets)
            for ifac,facet in enumerate(MP.facets):
                facet = list(trimSmall(array(facet)))
                temp = []
                for facetPoint in facet:
                    temp.append(facetPoint + point)
                tCell.facets[ifac] = temp
            self.facetsMathPrint(tCell,'v{}'.format(ipoint),False,'RandomColor[]')
            strOut = self.facetsMathToStr(strOut,tCell,'v{}'.format(ipoint),False,'RandomColor[]')
            print ';',; strOut+=';'
            showCommand += 'v{}'.format(ipoint)
            if ipoint < len(BZ.mesh)-1:
                showCommand += ','
        showCommand += ']'
        print ';',;strOut+=';'
        print showCommand ;strOut+=showCommand
        mathOut.write(strOut)
        mathOut.close()
    
    def mathPrintPoints(self,points):
        '''As spheres'''
        print 's = Graphics3D[{',
        for ipoint,point in enumerate(points):
            print 'Sphere[{' + '{},{},{}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(self.rpacking),
            if ipoint < len(points ) -1:
                print ',',
        print '}, Axes -> True, AxesLabel -> {"x", "y", "z"}]\n'

    def mathPrintPlanes(self,planes):
        print 'r=RegionPlot3D[',
        for iplane, pvec in enumerate(planes):
            print '{}x+{}y+{}z<={}&&'.format(pvec[0],pvec[1],pvec[2],dot(pvec,pvec)), #write plane equations for mathematica
        print ']\n'
        
        
#     def facetsMeshVCMathWrite(self,BZ,MP):
#         '''Put a mesh point VC at each mesh point'''
#         mathOut = open('mathOuts',w)
#         self.facetsMathPrint(BZ,'s','True','Red'); #draw supecell voronoi cell
#         showCommand = ';Show[s,'  
#         print ';',
#         for ipoint,point in enumerate(BZ.mesh):
#             tCell = cell()
#             tCell.facets = [[]]*len(MP.facets)
#             for ifac,facet in enumerate(MP.facets):
#                 facet = list(trimSmall(array(facet)))
#                 temp = []
#                 for facetPoint in facet:
#                     temp.append(facetPoint + point)
#                 tCell.facets[ifac] = temp
# #             self.facetsMathPrint(tCell,'v{}'.format(ipoint),False,\
# #                                  'discreteColors[{}][[{}]]'.format(len(BZ.mesh),ipoint));print ';',
#             self.facetsMathPrint(tCell,'v{}'.format(ipoint),False,'RandomColor[]');print ';',
# 
#             showCommand += 'v{}'.format(ipoint)
#             if ipoint < len(BZ.mesh)-1:
#                 showCommand += ','
# #                 print ',', 
#         showCommand += ']'
#         print ';',
#         print showCommand 
        
        
        
        
