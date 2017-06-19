'''
1. Get IBZ
2. Fill with an FCC mesh, with the desired number of points
3. Relax 100 steps with steepest descent

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

from symmetry import get_lattice_pointGroup, get_spaceGroup #these have vectors as ROWS

from kmeshroutines import (svmesh,svmeshNoCheck,svmesh1freedir, lattice_vecs, lattice, surfvol,
    orthdef, trimSmall, cosvecs,
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double,
    checksymmetry, nonDegen, MT2mesh, matchDirection,intoVoronoi,intoCell,
    reverseStructured,isInVoronoi,areParallel, among, addVec)
def timestamp():
    return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())

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

def directFromCart(Lvs,cartvec):
    '''Assumes lattice vectors are in **columns** in LVs.
    Then D = inv(transpose(LV)) * C'''
    return dot(inv((Lvs)), cartvec)
    
def cartFromDirect(self,Lvs,dirvec): 
    '''Assumes lattice vectors are in **columns** in LVs.
    Then C = transpose(LV) * D, because Cx = L1xD1 + L2xD2 + L3xD3, etc, 
    which is top row of transpose(LV) dot D, or the first matrix multiplication 
    operation'''
    return dot(Lvs, transpose(dirvec))

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
        [x2,y2,z2]] }  (ro^2,r1^2,r2^2) has a solution
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
    xy = planar3dTo2d(facet,eps)
    angles = []
    for i, vec in enumerate(facet):
        angle = arctan2(xy[i,1],xy[i,0])
        if angle < 0-eps: angle += 2*pi
        angles.append(angle)
    return [point for (angle,point) in sorted(zip(angles,facet),key = lambda x: x[0])] 

def planar3dTo2d(points,eps):
    '''Finds x,y coordinates of points in a plane'''
    coords = zeros((len(points),2))
    uvec = plane3pts(points,eps)[0] #plane normal
    rcenter = sum(points)/len(points)
    xunitv =  (points[0] - rcenter)/norm(points[0] - rcenter)
    crss = cross(xunitv, uvec)
    yunitv = crss/norm(crss)
    for i, vec in enumerate(points):
        vc = vec - rcenter
        coords[i,0] = dot(vc,xunitv); coords[i,1] = dot(vc,yunitv)
    return coords
    

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

def magGroup(arr,eps):
    '''Return a list of lists of indices, grouped by magnitude'''
    indsList = []
    magsList = arr['mag'].tolist()
    magsList.append(1000)
    lastBorder = 0
    for i,mag in enumerate(magsList):
        if i>0:
            if abs(magsList[i]-magsList[i-1]) > eps: #we have a border
                indsList.append(range(lastBorder,i))
                lastBorder = i
    return indsList


def newBounds(boundVecs,bndsLabels,grp,cell,type,eps):
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
        for pair in pairs:
            planes3 = [boundVecs[ig]['vec']]
            if not ig in pair:# or ig in keepLabels):
                planes3.append(boundVecs[pair[0]]['vec'])
                planes3.append(boundVecs[pair[1]]['vec'])
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
#     if mod(len(cell.bounds[0]),2)!=0:
#         sys.exit('Stop. Error in newBounds. BZ must have even number of planes')
    if len(cell.fpoints) >= 3 and type == 'BZ':
            checkNext = not areEqual(convexH(newVerts).volume,cell.volume,eps/4.0) 
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

def getVorCell(boundPlanesVecs,cell,type,eps):
    '''Boundaries and vertices of Voronoi cell'''   
    indsList = magGroup(boundPlanesVecs,eps)
    checkNext = True
    boundsLabels = []
#     for i in indsList[0]:  #add the first group
#         vec = boundPlanesVecs[i]['vec']; mag = norm(vec)
#         cell.bounds[0].append(vec/mag); cell.bounds[1].append(mag)
    for igroup, group in enumerate(indsList):
        #get any intersections between these planes
#         cell.fpoints = getInterscPoints(cell.bounds,eps)
        checkNext,boundsLabels,cell = newBounds(boundPlanesVecs,boundsLabels,group,cell,type,eps)
        if type == 'BZ' and not checkNext: break
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
        self.fpoints = [] #points as a set (but a list)
        self.volume = None
        self.center = None #body center of cell
        self.mesh = None #centers of each voronoi cell, or kpoints
        self.weights = None
    
class dynamicPack(): 
    ''''''
    from comMethods import readfile,writefile,trimSmall,areEqual,directFromCart,cartFromDirect
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
    from conjGradMin2 import (fmin_cg,minimize_cg,line_search_wolfe1,scalar_search_wolfe1)
        
    def pack(self,A,B,totatoms,aTypes,postype,aPos,targetNmesh,meshtype,path,method):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        self.B = B
#         print '\nB (Recip lattice vectors as columns',B
#         print 'method',method
        self.initFactor = 1.0 #assign 90% of the point at the beginning based on an FCC or BCC mesh
        self.power = 6.0
        self.wallfactor = 1.0  #probably needs to be bigger than interfactor by about the average number of nearest neighbors
        self.wallClose = 0.5 #0.5 #to allow initial points closer to the wall set to less than 1. 
        self.wallOffset = 0.75 #back off wall forces and energies by a distance that is a fraction of dw. 
        self.interfactor = 1.0        
        self.nTarget = int(self.initFactor*targetNmesh)
        self.path = path
        self.method = method
            #0: exact: use vertices of mesh voronoi cell that are closest/farthest 
            #         from the IBZ center origin to check if the point's volume is cut. 
            #         Cut the VC to determine the volume contribution      
            #0.5 approx  If cut volume is less than 50%, distribute weight to neighbors of equivalent points
            #         If point is outside of first BZ, then translate by G vector to get it inside.  Points inside but not in IBZ
            #        use point symmetries to get in.  This applies to kpoints in IBZ but near corners as well. 
        vol = abs(det(B))
        self.ravg = (vol/targetNmesh)**(1/3.0) #distance if mesh were cubic. 
        self.df = self.ravg #inter-point force scale distance
#         self.dw = self.df/2.0 #wall force scale distance
        self.dw = self.df/2 #wall force scale distance
        self.shift =  array([1,1,1])/8.0 #array([1/10,0,0])
        eps = self.ravg/2000
#         self.symops = get_lattice_pointGroup(transpose(self.B),eps)
#         self.nops = len(self.symops)
        [symopsList, fracsList] = get_spaceGroup(transpose(A),aTypes,transpose(aPos),1e-8,postype.lower()[0] == 'd')
        self.nops = len(symopsList)
        self.symops = zeros((3,3,self.nops),dtype = float)
        for iop in range(len(symopsList)):
            self.symops[:,:,iop] = array(symopsList[iop])
        
#         get_spaceGroup(par_lat,atomType,bas_vecs,eps=1E-10,lattcoords = False):
#           ...
#         return(sg_ops,sg_fracts)
#     Args:
#         par_lat (array-like): A 2D integere array that contains the parent lattice vectors
#         atomType (list of int): Integer array representing the type of each basis atom
#         bas_vecs (array-like): A 2D integere array that contains the basis vectors for the cell
#         eps (float, optional): Finite precisions tolerance
# 
#         lattcoords (bool, optional): True if vectors are in lattice coordinates 
#           rather than cartesian
# 
#     Returns:
#         sg_ops (array-like): The rotation and mirror operations of the space group.
#         sg_fracts (array-like): The translation operations of the space group.
        
        print 'Number of desired points:', targetNmesh
        print 'Symmetry operations:', self.nops

        BZ = cell() #instance
        BZ.volume = vol
        braggVecs = getBraggVecs(self.B)
        BZ = getVorCell(braggVecs,BZ,'BZ',eps)
#         print 'Vornonoi cell'; self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n'
        self.vorCell = BZ
#         self.facetsMathPrint(BZ,'p',True,'Red') 
        IBZ = self.getIBZ(BZ,eps) #now irreducible BZ
#         self.facetsMathPrint(IBZ,'p',True,'Red');print ';Show[p]\n' 
        IBZ = self.meshInitCubic(IBZ,meshtype,eps)
        if 0 < len(IBZ.mesh) <= 1000:
            OK = True
            self.dynamic(IBZ,eps)
            IBZ = self.weightPoints(IBZ,eps)
            self.writeKpoints(IBZ)
            return OK
        else: 
            OK = False
            return OK
            
#         sys.exit('stop')
        return
    
    def weightPoints(self,IBZ,eps):
        '''Find the volume of the Voronoi cell around each point, and use it to weight the point.
        Search a sphere of radius a few df for neighbors.  Use the half vectors to these points 
        and the vectors to the walls to define the bounding planes.
        
        Vectors need to be taken from each mesh point as the origin'''
                
        for ip,point in enumerate(IBZ.mesh):
            pointCell = cell()
            neighs,neighLbls = self.getNeighbors(point,IBZ,eps)
            boundVecs = zeros(len(neighs)+ len(IBZ.bounds[0]),dtype = [('vec', '3float'),('mag', 'float')]) 
            for iw, u in enumerate(IBZ.bounds[0]):
                ro = self.bounds[1][iw]
                d = ro-dot(point,u) 
                vec = d*u
                boundVecs[iw]['vec'] = vec  #vector from point to wall
                boundVecs[iw]['mag'] = norm(vec)
            for j, jpoint in enumerate(neighs):
                vec = (jpoint - point)/2
                boundVecs[j+len(IBZ.bounds[0])]['vec'] = vec
                boundVecs[j+len(IBZ.bounds[0])]['mag'] = norm(vec)
            boundVecs.sort(order = 'mag')
            pointCell = getVorCell(boundVecs,pointCell,'point',eps)
#             print 'volume',pointCell.volume, 'number of these to fill IBZ', IBZ.volume/pointCell.volume
#             IBZ.weights.append(pointCell.volume/self.ravg**3) 
            IBZ.weights.append(pointCell.volume) 
#             print 'Point vor cell'; self.facetsMathPrint(pointCell,'p',True,'Red');print ';Show[p]\n'

        wtot = sum(IBZ.weights)
        print 'Total volume of point Vor cells',wtot,'vs IBZ volume', IBZ.volume
        if not areEqual(wtot, IBZ.volume, eps):
            sys.exit('Stop: point Voronoi cells do not sum to the IBZ volume.')
        IBZ.weights /= self.ravg**3  #to scale them order(1).  
        return IBZ

    def meshInitCubic(self,IBZ,type,eps):
        '''Add a cubic mesh to the interior, . If any 2 or 3 of the facet planes are 
        orthogonal, align the cubic mesh with their normals.       
        Remove any points within self.rw from any wall        '''
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
                vec = pairs[-1][i]
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
        #Find the extremes in each cubLV direction:
        intMaxs = [] #factors of aKcubConv
        intMins = []
        for i in range(3):
            projs = []
            for point in IBZ.fpoints:
                shifted = point + self.shift
                projs.append(dot(cubicLVs[:,i],shifted)/aKcubConv**2)
            intMaxs.append(int(ceil(max(projs)))+2) #optimize: Why is +2 required with shift of 1/2,1/2,1/2 on cubic?
            intMins.append(int(floor(min(projs)))-1)#optimize: Is -1 required? 
#         print 'Maxes',intMaxs
#         print 'Mins',intMins       
        #Create the cubic mesh inside the irreducible BZ
        IBZ.mesh = []
        IBZ.weights = []
        weightsInside = 0
        nInside = 0         
        ik = 0 
        #begin MP facets printing
#         self.facetsMathPrint(IBZ,'s','True','Red'); print ';', #draw supecell voronoi cell before loop
        showCommand = 'Show[s,' 
        #done with start of MP facets printing        
        for i in range(intMins[0],intMaxs[0]):
            for j in range(intMins[1],intMaxs[1]):
                for k in range(intMins[2],intMaxs[2]):
                    lvec = i*cubicLVs[:,0]+j*cubicLVs[:,1]+k*cubicLVs[:,2]
                    for iS, site in enumerate(sites):
                        ik+=1
                        kpoint = lvec + aKcubConv*self.shift + site
#                         print 'kp',kpoint
                        if isInside(kpoint,IBZ.bounds,self.dw*self.wallClose):  #Can't be closer than self.dw*self.wallClose to a wall
#                         if i + j + k ==1:
#                             'pause'
#                         if not isOutside(kpoint,IBZ.bounds,eps):  #Can't be closer than self.dw*self.wallClose to a wall
                            nInside += 1
                            IBZ.mesh.append(kpoint)
#                             IBZ.weights.append(1.0) 
        print 'Points inside', nInside
#         if nInside == 0: sys.exit('No points are inside the IBZ.  Increase the number of target points')
#         print 'BZ volume:', det(self.B),'\n'
#         self.facetsMeshMathPrint(IBZ); print ';Show[p,q]\n'
        return IBZ

    def dynamic(self,IBZ,eps):
        ''' '''
        self.IBZ = IBZ
        self.points = IBZ.mesh #initial mesh
        self.bounds = IBZ.bounds 
        self.facets = IBZ.facets
        self.eps = eps
        self.relax()

    def relax(self):
        '''Minimization of the potential energy.
        The energy must be a function of 1-D inputs, so we flatten the points into components '''
        initialPos = array(self.points).flatten()
        epsilon = self.ravg/100
#         out = self.fmin_cg(initialPos,self.eps)
        self.minNewtSteepest(initialPos,self.eps)    
        return

    def enerGrad(self,comps):
        '''Returns the total energy, gradient (-forces) and (optional) Hessian matrix, 
        using comps, which are the flattened components of the current positions  '''
#         print 'oldindvecs',self.oldindVecs
        self.forces = zeros((len(self.points),3))
        self.wallForce = zeros(len(self.IBZ.facets))
        self.wallPress = zeros(len(self.IBZ.facets))
        vecs = comps.reshape((len(self.points),3))
        p = self.power
        wallfact = self.wallfactor
        interfact = self.interfactor
#         self.getHessian(wallfact,interfact,p)
        etot = 0
        for i,ri in enumerate(vecs):
#             print 'EG{}'.format(i),
#             print 'point',i,ri
            #wall forces
            for iw, u in enumerate(self.bounds[0]):
                ro = self.bounds[1][iw]
                d = ro-dot(ri,u)+ self.wallOffset*self.dw #distance from plane to ri offset factor allows points to move closer to walls. 
#                 if d<0:
#                     'pause'
# #                     print '** point {} crossed boundary'.format(i), iw, u,ro
# #                     d = -d/10 #Have crossed boundary. Increase force by shortening effective d. 
                if d<0:
                    print 'ri,ro,u, dot(ri,u),d'
                    print ri,ro,u, dot(ri,u), d 
                    sys.exit('Error. Point {} in enerGrad is not in the IBZ.'.format(iw))
                fmag = wallfact*(d/self.dw)**(-p)  #dimensionless
                etot += wallfact*self.dw/abs(-p+1)*(d/self.dw)**(-p+1)#Units of length. Both F and E can't be dimensionless unless we make the positions dimensionless.
#                 print '\t wall',iw,d, u,ro
#                 print '\t\tf',-u*fmag,'\te',wallfact*self.dw/abs(-p+1)*(d/self.dw)**(-p+1)
                self.forces[i] += -u*fmag
                self.wallForce[iw] += fmag #since forces are normal to plane, we sum the magnitudes
#             print '\tfwtot',self.forces[i]
#                 print 'Wall',iw,d,'force',-u*fmag,fmag
#            inter-point forces
            for j, rj in enumerate(vecs):
                if i!=j:
                    d = norm(ri-rj)
#                     print 'Inter d,f', d,interfact*(d/self.df)**(-p)*(ri-rj)/d
                    self.forces[i] += interfact*(d/self.df)**(-p)*(ri-rj)/d
                    if j>i: #don't overcount
                        etot += interfact*self.df/abs(-p+1)*(d/self.df)**(-p+1)
        for i,fac in enumerate(self.facets):
            area = convexH(planar3dTo2d(fac,self.eps)).volume  # for 2d problems, the "volume" returned is the area, and the "area" is the perimeter
            self.wallPress[i] = self.wallForce[i]/area
#         print '\tetot,norm',etot,norm(-self.forces.flatten())
        return etot, -self.forces.flatten() #gradient is opposite the forces.
        
    def writeKpoints(self,cell):
        nk = len(cell.mesh)
        totw = sum(cell.weights)
        lines = []
        lines.append('Packing of IBZ (Bret Hess, BYU). Total weights: {:12.8f} \n'.format(totw))#(vs 1.0 per general point without symmetry\n'.format(totw))
        lines.append('{}\n'.format(nk))
        lines.append('Reciprocal\n') #direct coordinates!...vasp doesn't read cartesian kpoints right
        for ik,kpoint in enumerate(cell.mesh):
            #direct coordinates!...vasp doesn't read cartesian kpoints right
            kpointDir = directFromCart(self.B,kpoint)
            lines.append('{:15.12f}  {:15.12f}  {:15.12f}  {:15.12f}\n'\
                         .format(kpointDir[0],kpointDir[1],kpointDir[2],cell.weights[ik]))
        self.writefile(lines,'KPOINTS')         
       
    def intoIBZ(self,kpoint,IBZ,eps):   
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
                    
    def getNeighbors(self,kpoint,IBZ,eps):
        '''Search a sphere around the kpoint and collect neighbors.
        Then if outside the IBZ, move point into IBZ by symmetry.  Search another sphere.
        Return the neighbors
        '''
        neighR = 8.0*self.rpacking
        neighs,neighLbls = self.searchSphere(kpoint,IBZ,neighR,eps)
        return neighs,neighLbls 
    
    def searchSphere(self,kpoint,IBZ,R,eps):
        neighs = []
        neighLbls = []
        for ip,meshPt in enumerate(IBZ.mesh):
            if not allclose(kpoint,meshPt,eps):
                if norm(meshPt-kpoint) <= R:
                    neighs.append(meshPt)
                    neighLbls.append(ip)
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
        for ifac, facet in enumerate(cell.facets):
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
                        interscP = intsPlLinSeg(u,ro,pointi,pointj,eps) #         
                        newfacet.append(interscP)
                        bordersFacet = addVec(interscP,bordersFacet,eps)

                if len(newfacet) >= 3:
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
#         print '\n\nReducing Brillouin zone by symmetry'
#         self.facetsMathPrint(BZ,'p',True,'Red');print ';Show[p]\n' 
#          
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
                        #the plane to cut is the plane of O and axis, so take normal perpendicular to vector O.                   )
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
                    
                    
#                     BZ.volume = convexH(BZ.fpoints).volume
#                     self.IBZvolCut = det(self.B)/BZ.volume
#                     getBoundsFacets(BZ,eps)
#                     BZ.fpoints = flatVecsList(BZ.facets,eps)
#                     BZ.center = sum(BZ.fpoints)/len(BZ.fpoints)
#                     print 'Vol BZ / Vol IBZ', self.IBZvolCut                    
                    
                    
                    
                    
                      
            elif areEqual(det(op),-1.0,eps):
                inversion = True
        if inversion and self.nops==2: #apply last of all.  For now I think inversion acts on the IBZ only if it is alone with the id
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
        if not areEqual(self.IBZvolCut,self.nops,eps):
            sys.exit('Volume not reduced by factor equal to the number of symmetry operations')
        return BZ

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
        '''Writes facets and spheres for Mathematica drawing'''
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
        self.facetsMathPrint(BZ,'s','True','Red'); #draw supercell voronoi cell
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
        
    def plotPos(self,arr,npoints,step):
        self.plot2dPts(arr,npoints,timestamp()+step,'x','y') 
        self.plot2dPts(arr,npoints,timestamp()+step,'x','z')  
        
    def plot2dPts(self,arr,npoints,tag,ax0,ax1):
#         fig = figure()
        i0 = ('x','y','z').index(ax0)
        i1 = ('x','y','z').index(ax1)
        print 'arr0',arr[:,i0]
        print 'arr1',arr[:,i1]
        figure()
        scatter(arr[:,i0],arr[:,i1])
        xlabel('{}'.format(ax0))
        ylabel('{}'.format(ax1))
        name = '{}_{}-{}'.format(tag,ax0,ax1)
        title(name)
        savefig(self.path+'/'+ name+'.pdf');
        os.system('cp {} {}'.format(self.path+'/'+ name+'.pdf',self.path+ '/latest{}-{}'.format(ax0,ax1)+'.pdf' ))
        close()
        
    def minNewtSteepest(self,x0,eps):
        ''' See intro to P. Pulay, Chem. Phys. Lett. 73, 393 (1980) for use of Hessian.  
        But a steepest descent works better on first tests (see en.wikipedia.org/wiki/Gradient_descent). 
        We move in the direction of -gradient, starting with a default step, and decreasing it by 2 until
        a lower energy and lower norm(gradient) is found.   The default step is increased if the previous step succeeded without adjustment.
        If the step must be lowered to find a lower energy and lower norm(gradient), the default is lowered. 
        
        '''
        self.IBZ.mesh = self.points
#         self.facetsMeshMathPrint(self.IBZ); print ';Show[p,q]\n'
#         print 'Start Minimization';

        itermax = 100
        gnormTol = 0.01
        minstep = 0.0001
        xold= x0
        fold,gold = self.enerGrad(xold)
        fnew = fold
        gnormold = norm(gold)
        fstart = fold; gstart = gold; gnormstart = gnormold
        method = 'steepest'
#         if method =='newtRaph': H = self.getHessian();Hinv = inv(H)    
#         H = self.getHessian();Hinv = inv(H)    
        
#         xolder =  xold + dot(Hinv,gold) #go backwards
        xolder =  xold + 0.01*gold #go backwards
#         golder = dot(H,(xolder-xold))
        folder = dot(gold,(xolder-xold))
        print 'energy_0',fold, 'gnorm',gnormold #,'grad', gold
        fstart = fold; gnormstart = gnormold
        iIter = 0
        step = 1.0
        atMinStep = False
        while iIter < itermax and gnormold > gnormTol and not atMinStep:
            print iIter, #progress bar
            method = 'steepest'
#             if method == 'steepest':
#                 step = max(step,dot(xold-xolder,gold-golder)/norm(gold-golder))  #barzilai-Borwein method                   
            lower = False
            while not lower:
                if step < minstep:
                    print 'minimum step reached: {}'.format(step) 
                    atMinStep = True
                    break
#                     if method == 'newtRaph':
#                         self.IBZ.mesh = self.points; self.facetsMeshMathPrint(self.IBZ); print ';Show[p,q]\n'
#                         sys.exit('Step {} is still too small after "steepest" atempt: stop'.format(step))
#                     #try steepest descent
#                     method = 'newtRaph'
# #                     if method == 'steepest':
# #                         step = dot(xold-xolder,gold-golder)/norm(gold-golder)  #barzilai-Borwein method
# #                         if step < minstep:
# #                             step = 1.0
#                     print 'try newtRaph'              
#                 print 'step',step
#                 if method == 'newtRaph':
#                     xnew = xold - step*dot(Hinv,gold)
                if method == 'steepest':
                    xnew = xold - step*gold
                self.points = xnew.reshape((len(self.points),3))
                inside = True
                for point in self.points:
                    if not isInside(point,self.bounds,eps):
                        print 'point is outside IBZ...reduce step size:',point
                        inside = False
                        break                   
                if inside:
                    fnew,gnew = self.enerGrad(xnew)
                    gnormnew = norm(gnew)
#                     print '\nenergy',iIter, fnew,'step',step,'gnorm',gnormnew #, 'grad', gnew
                    if fnew<fold and gnormnew < gnormold:
                        lower = True
                step /= 2
            step *= 4
            xolder = xold
            folder = folder
            golder = gold
            xold = xnew
            fold = fnew
            gold = gnew
            gnormold = gnormnew
            iIter += 1
#             if method =='newtRaph': H = self.getHessian();Hinv = inv(H)                        
        self.IBZ.mesh = self.points; 
        print;self.facetsMeshMathPrint(self.IBZ); print ';Show[p,q]\n'
        print 'For {} points in IBZ and {} steps'.format(len(self.points),iIter)
        print '\tStarting energy',fstart, 'gnorm',gnormstart
        print '\tEnding energy',fnew,'gnorm',gnormnew, 'step',step,#, 'grad', gnew
        if gnormnew <= gnormTol:
            print '\nSuccess after {} iterations'.format(iIter)
        elif iIter == itermax:
            print '\nExceeded maximum number of iterations ({}), while gnorm {} is greater than the tolerance {}'.format(itermax,gnormnew,gnormTol)
        if not (fnew < fstart and gnormnew < gnormstart):
            sys.exit('Did not find a lower energy and force norm: stop')
#         print; print
        
        
    def getHessian(self):
        '''
        For a wall force in the form: f = (d/dw)^(-p), 
        D[en[d[x]], x, y] = -p ux uy (d/w)^(-1 - p)/w
       
        For an interparticle force in the form: f = (d/df)^(-p), the pair energy is 
        of the form 
        *Same type* component (x or y or z) on *different particles*:  Here x, y z represent any three components,
        reordered here:
          D[en[d[x1, x2, y1, y2]], x1, x2] =  
            -(d/df)^-p/d^3 *(p(x1 - x2)^2 - y1^2 + 2 y1 y2 - y2^2 - z1^2 + 2 z1 z2 - z2^2)
            For example, for particles a, b, c,  d^2E/dx_a dx_b
            This applies to the one pair a-b
        
        Identical component on the *same* particle: D[en[d[x1, x2, y1, y2]], x1, x1] = -(above result)
            For example, for particles a, b, c,  d^2E/(dx_a)^2
            This applies to all pairs that include particle a
        
        *Different type* components (x vs y vs z) on same or different particles (x vs y vs z): 
          D[en[d[x1, x2, y1, y2]], x1, y1] = (d/df)^-p/d^3 * (1 + p)(x1 - x2)(y1 - y2)
            For example, for particles a, b, c,  d^2E/dx_a dy_b or d^2E/dx_a dy_a 
            If comps are on the *same particle*, this applies to all pairs that include particle a 
            If comps are on *different particles*, this applies to the one pair a-b.  But use y2-1'''
        p = self.power
        wallfact = self.wallfactor
        interfact = self.interfactor
        points = self.points
        N = len(points)
        self.hessian = zeros((3*N,3*N),dtype = float)
        df = self.df
        dw = self.dw 
        for ipoint in range(N): #differentiation particle
            for icomp in range(3): #differentiation comp
                for jpoint in range(ipoint,N): #2nd differentiation particle (not pair!)
                    for jcomp in range(3): #differentiation comp
#                         print '\ndiff',ipoint,jpoint,'comps',icomp,jcomp
                        if ipoint == jpoint:
                            #wall forces: only same-particle contributes, but could be different comps
                            for iw, u in enumerate(self.bounds[0]):
                                ro = self.bounds[1][iw]
                                d = ro-dot(points[ipoint],u) #distance from plane to ri
                                
                                if d<0:
    #                                 d = -d/10 #Have crossed boundary. Increase force by shortening effective d. 
                                    print 'ri,ro,u, dot(ri,u),d'
                                    print ri,ro,u, dot(ri,u), d
                                    sys.exit('Error. Point {} in getHessian is not in the IBZ.'.format(iw))
#                                 print 'Wall',wallfact*(p*u[icomp]*u[jcomp]*(d/dw)**(-1 - p))/dw, u
                                self.hessian[ipoint*3+icomp,ipoint*3+jcomp] += wallfact*(p*u[icomp]*u[jcomp]*(d/dw)**(-1 - p))/dw
                            #differentiation coordinates on same particles
                            for mpoint in range(ipoint+1,N):  #forming all energy pairs that contain ipoint
                                comps = [0,1,2]
                                d = norm(points[ipoint]-points[mpoint])
                                x1 = points[ipoint][icomp]
                                x2 = points[mpoint][icomp]
                                comps.remove(icomp)
                                if icomp == jcomp:
                                    y1 = points[ipoint][comps[0]]
                                    y2 = points[mpoint][comps[0]]                 
                                    z1 = points[ipoint][comps[1]]
                                    z2 = points[mpoint][comps[1]]
#                                     print 'pair',ipoint,mpoint,interfact*((d/df)**-p)/d**3 *(p*(x1 - x2)**2 - y1**2 + 2*y1*y2 - y2**2 - z1**2 + 2*z1*z2 - z2**2)
                                    self.hessian[ipoint*3+icomp,ipoint*3+icomp] += interfact*((d/df)**-p)/d**3 *(p*(x1 - x2)**2 - y1**2 + 2*y1*y2 - y2**2 - z1**2 + 2*z1*z2 - z2**2)
                                else:
                                    y1 = points[ipoint][jcomp]
                                    y2 = points[mpoint][jcomp]
                                    comps.remove(jcomp)
                                    z1 = points[ipoint][comps[0]]
                                    z2 = points[mpoint][comps[0]]
#                                     print 'other mpoint',mpoint,interfact*((d/df)**-p)/d**3 * (1 + p)*(x1 - x2)*(y1 - y2)      
                                    self.hessian[ipoint*3+icomp,ipoint*3+jcomp] += interfact*((d/df)**-p)/d**3 * (1 + p)*(x1 - x2)*(y1 - y2)                
                        else: #different particles, which must be ipoint and jpoint
                            comps = [0,1,2]
                            d = norm(points[ipoint]-points[jpoint])
                            x1 = points[ipoint][icomp]
                            x2 = points[jpoint][icomp]
                            comps.remove(icomp)
                            if icomp == jcomp: #same comp type
                                y1 = points[ipoint][comps[0]]
                                y2 = points[jpoint][comps[0]]                 
                                z1 = points[ipoint][comps[1]]
                                z2 = points[jpoint][comps[1]]
#                                 print 'i-j interpoint', -interfact*((d/df)**-p)/d**3 *(p*(x1 - x2)**2 - y1**2 + 2*y1*y2 - y2**2 - z1**2 + 2*z1*z2 - z2**2)
                                self.hessian[ipoint*3+icomp,jpoint*3+jcomp] += -interfact*((d/df)**-p)/d**3 *(p*(x1 - x2)**2 - y1**2 + 2*y1*y2 - y2**2 - z1**2 + 2*z1*z2 - z2**2)
                            else: #different comp type
                                y1 = points[ipoint][jcomp]
                                y2 = points[jpoint][jcomp]
                                comps.remove(jcomp)
                                z1 = points[ipoint][comps[0]]
                                z2 = points[jpoint][comps[0]]
#                                 print 'i-j interpoint', interfact*((d/df)**-p)/d**3 * (1 + p)*(x1 - x2)*(y2 - y1) 
                                self.hessian[ipoint*3+icomp,jpoint*3+jcomp] += interfact*((d/df)**-p)/d**3 * (1 + p)*(x1 - x2)*(y2 - y1) 
        for i in range(3*N):
            for j in range(i,3*N):
                self.hessian[j,i] = self.hessian[i,j]
                                                        
        return self.hessian   
                    
                    
                    
