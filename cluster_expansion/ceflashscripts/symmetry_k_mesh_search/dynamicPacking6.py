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
                   trace,where,real,allclose,sign,pi,imag,identity,argsort,
                   cos,sin,log)
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
from scipy.spatial import Voronoi
from __builtin__ import True
from _symtable import CELL
# import datetime
# from _ast import operator
# from pip._vendor.html5lib.constants import rcdataElements
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts')
sys.path.append('/fslhome/bch/graphener/graphener')

from symmetry import get_lattice_pointGroup, get_spaceGroup #these have vectors as ROWS

def readfile(filepath):
    file1 = open(filepath,'r')
    lines = file1.readlines()
    file1.close()
    return lines

def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

def areEqual(x,y,eps):
    return abs(x-y)<eps

def isinteger(x,eps):
    return areEqual(abs(rint(x)-x),0,eps)

def icycle(i,change): #for cycling indices 0,1,2
    i = i+change
    if i>2:
        i=0
    if i<0:
        i=2
    return i

def trimSmall(list_mat):
    low_values_indices = abs(list_mat) < 1.0e-5
    list_mat[low_values_indices] = 0.0
    return list_mat

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
    magsList.append(1000) #need to end the list with a final number to stop on. 
    lastBorder = 0
    for i,mag in enumerate(magsList):
        if i>0:
            if abs(magsList[i]-magsList[i-1]) > eps: #we have a border
                indsList.append(range(lastBorder,i))
                lastBorder = i
    return indsList


def newBounds(boundVecs,bndsLabels,grp,cell,type,eps):
    '''Find intersection points:  
    If any of these intersections are inside the current boundaries, then they 
    become new vertices, and the planes are added to boundaries. 
    
    If the new volume is the same as det(B), then we are done
    '''
    keepLabels = []
    checkNext = True
    allPlanes = [cell.bounds[0][i]*cell.bounds[1][i] for i in range(len(cell.bounds[0]))]
    allVerts = deepcopy(cell.fpoints)
    for ig  in grp:
        bndsLabels.append(ig)
    pairs = list(combinations(bndsLabels,2))
    for ig in grp:
        for pair in pairs:
            planes3 = [boundVecs[ig]['vec']]
            if not ig in pair:# or ig in keepLabels):
                planes3.append(boundVecs[pair[0]]['vec'])
                planes3.append(boundVecs[pair[1]]['vec'])
                intersPt = threePlaneIntersect(planes3)   
                if not intersPt is None:
                    if not isOutside(intersPt,cell.bounds,eps)\
                        and not among(intersPt,allVerts,eps):
                            addVec(intersPt,allVerts,eps)
                            addVec(planes3[0],allPlanes,eps)
                            addVec(planes3[1],allPlanes,eps)
                            addVec(planes3[2],allPlanes,eps)                                
                                     
    if len(allVerts)>0:
        #keep only the vertices that can be reached without crossing any plane
        newVerts = []
        for vert in allVerts:
            for plane in allPlanes:
                if dot(vert,plane)>dot(plane,plane)+eps:
                    break
            else: 
                if norm(vert)>30:
                    'pause'
                newVerts.append(vert)
        #Start new to remove planes that don't host a vertex.
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
#     self.mathPrintPlanes(allPlanes)
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
    cell.facets = []
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
        if len(facetvecs)> 3:          
            cell.facets.append(orderAngle(facetvecs,eps))
        elif len(facetvecs)==3:
            cell.facets.append(facetvecs)
#             sys.exit('stop. Plane {} has less than three facetvecs')
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
#             print points[i],i,'on axis or mirror'
            break #no rotation effect
        otherLabels = range(len(points))
        otherLabels.pop(i)
        for j in otherLabels:
            if allclose(rpoint,points[j],atol=eps):
#                 print points[i],i,'maps to',j,rpoint
                return True,rpoint,points[j] 
#         else:
#             print points[i],i,'no map'
    return False,'',''

# def makesAllDups(op,cell,eps):
#     '''Applies symmetry operator to all facet points. If all facet points are 
#     moved on top of other facet points, then return true'''
#     dups = zeros(len(cell.fpoints),dtype = bool)
#     points = cell.fpoints
#     print len(points),'points'
#     for i in range(len(points)):
#         rpoint = dot(op,points[i])
#         if allclose(points[i],rpoint,atol=eps):
# #             print points[i],i,'on axis or mirror'
#             dups[i] = True
#             break #no rotation effect...point is on axis or mirror
#         otherLabels = range(len(points))
#         otherLabels.pop(i)
#         for j in otherLabels:
#             if allclose(rpoint,points[j],atol=eps):
# #                 print points[i],i,'maps to',j,rpoint
#                 dups[i] = True
#                 break
# #         else:
# #             print points[i],i,'no map'
#     return all(dups)

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

def intsPlLinSeg(u,ro,r1,r2,eps):
    '''Intersection between a plane and a line.
    A plane is given by dot(r,u) = ro.  
    A line segment between r1 and r2 is given by vectors r = r1+t(r2-r1), for t in (0,1) 
    Combining these:  r1u = dot(r1,u).  r2u = dot(r2,u).  Then t = (ro-r1u)/(r2u-r1u).
    So if t is in (0,1), then we have an intersection'''
    r1u = dot(r1,u)
    r2u = dot(r2,u)
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
    return cell

class cell():
    def __init__(self):
        self.bounds = [[],[]] #planes, written as normals and distances from origin [u's] , [ro's]
        self.facets = None #points arranged in facets
        self.fpoints = [] #points as a set (but a list)
        self.volume = None
        self.center = None #body center of cell
        self.mesh = None #centers of each voronoi cell, or kpoints
        self.weights = None
    
class dynamicPack(): 
    ''''''
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
    from conjGradMin2 import (fmin_cg,minimize_cg,line_search_wolfe1,scalar_search_wolfe1)
        
    def pack(self,A,B,totatoms,aTypes,postype,aPos,targetNmesh,meshtype,path,method):
        #1. nCoarse random points in the cell parallelpiped.  
#         nCoarseMax = 200
        self.B = B
#         print '\nB (Recip lattice vectors as columns',B
#         print 'method',method
        vol = abs(det(B))
        self.ravg = (vol/targetNmesh)**(1/3.0) #distance if mesh were cubic. 
        
        self.power = 6.0
        self.wallPower = 6.0
        self.wallfactor = 1.0  #probably needs to be bigger than interfactor by about the average number of nearest neighbors
        self.wallClose = 0.1 #0.5 #to allow initial points closer to the wall set to less than 1. 
        self.wallOffset = 0.5 #back off wall forces and energies by a distance that is a fraction of dw. 
        self.interfactor = 1.0        
        self.initFactor = 1.0
        self.df = 1.00 * self.ravg #inter-point force scale distance
        self.dw = 0.5 * self.df #wall force scale distance
#        self.shift =  array([1,1,1])/8.0 #array([1/10,0,0])
        eps = self.ravg/300
        self.eps = eps
        self.searchInitFactor = 0.0
        self.initSrch = 'max'
#         self.initSrch = None
#         self.initSrch = 'target'
        self.nTarget = int(self.initFactor*targetNmesh)
        self.path = path
        self.method = method
        [symopsList, fracsList] = get_spaceGroup(transpose(A),aTypes,transpose(aPos),1e-3,postype.lower()[0] == 'd')
        self.nops = len(symopsList)
        self.symops = zeros((3,3,self.nops),dtype = float)
        for iop in range(len(symopsList)):
            self.symops[:,:,iop] = trimSmall(array(symopsList[iop]))       
        print 'Number of desired points in full BZ:', targetNmesh
        BZ = cell() #instance
        BZ.volume = vol
        braggVecs = getBraggVecs(self.B)
        BZ = getVorCell(braggVecs,BZ,'BZ',eps)
        self.facetsMathFile(BZ,'BZ') 
        self.IBZ = self.getIBZ(BZ,eps) #now irreducible BZ
        self.nTargetIBZ = int(rint(self.nTarget/float(self.nops)))
        self.facetsMathFile(self.IBZ,'IBZ') 
        self.meshInitCubic(meshtype,eps)
        if 0 < len(self.IBZ.mesh) <= 4000:
            OK = True
            self.dynamic(eps)
            self.weightPoints(eps)
            self.writeKpoints()
            self.writeSym()
            return OK,self.nops
        else: 
            OK = False
            return OK,self.nops
    
    def sortfpoints(self,cell):
        #sort with vertices closest to the origin first
        norms = []
        newfpoints = zeros(len(cell.fpoints))
        for point in cell.fpoints:
            norms.append(norm(point))
        order = argsort(norms)
#         for inew,iold in order:
#             newfpoints.append(cell.fpoints[iold])
#         cell.fpoints = newfpoints
        cell.fpoints = array(cell.fpoints)[order].tolist()
        return cell
            
    def writeSym(self):
        writefile(['nops: {}\n'.format(self.nops),\
                   'IBZvolCut: {}\n'.format(self.IBZvolCut),\
                   'IBZvol: {}\n'.format(self.IBZ.volume)],'sym.out')

    def weightPoints(self,eps):
        '''Find the volume of the Voronoi cell around each point, and use it to weight the point.
        Search a sphere of radius a few df for neighbors.  Use the half vectors to these points 
        and the vectors to the walls to define the bounding planes.
        
        Vectors are first taken from each mesh point as the origin, 
        then displaced to their real positions in the cell for possible display'''
        allMPfacets = []
        for ip,point in enumerate(self.IBZ.mesh):
            print ip,
            pointCell = cell()
            neighs,neighLbls = self.getNeighbors(point,self.IBZ,eps)
#             print 'neighLbls',neighLbls
            boundVecs = zeros(len(neighs)+ len(self.IBZ.bounds[0]),dtype = [('vec', '3float'),('mag', 'float')]) 
            for iw, u in enumerate(self.IBZ.bounds[0]):    
                ro = self.IBZ.bounds[1][iw]
                d = ro-dot(point,u) 
                vec = d*u
                boundVecs[iw]['vec'] = vec  #vector from point to wall
                boundVecs[iw]['mag'] = norm(vec)
#                 print 'wall',iw,u, vec, norm(vec)
            for j, jpoint in enumerate(neighs):
                vec = (jpoint - point)/2
                boundVecs[j+len(self.IBZ.bounds[0])]['vec'] = vec
                boundVecs[j+len(self.IBZ.bounds[0])]['mag'] = norm(vec)
#                 print 'neighs',j,jpoint, vec, norm(vec)
            boundVecs.sort(order = 'mag')
            pointCell = getVorCell(boundVecs,pointCell,'point',eps)
            self.IBZ.weights.append(pointCell.volume)
            
            #For completeness,could update pointCell.center and pointCell.fpoints.  For brevity, we don't do this. 

            allMPfacets.append(pointCell.facets)

        print
        self.facetsMeshVCMathFile(self.IBZ,allMPfacets)
        wtot = sum(self.IBZ.weights)
        stdev = std(self.IBZ.weights)
        meanV = mean(self.IBZ.weights)
        volCheck = 0.1
        volErr = wtot - self.IBZ.volume        
        volErrRel = volErr/self.IBZ.volume
        
        print 'Total volume of point Vor cells',wtot,'vs IBZ volume', self.IBZ.volume
        print 'Relative volume error', volErrRel,'Abs volume error', volErr, 'Std dev/mean',stdev/meanV
        if not areEqual(wtot, self.IBZ.volume, volCheck*self.IBZ.volume):
#             print 'Total volume of point Vor cells',wtot,'vs IBZ volume', self.IBZ.volume
            sys.exit('Stop: point Voronoi cells do not sum to the IBZ volume.')
        else:
            print 'Point Voronoi cells volumes sum OK to within factor of {} of IBZ volume OK'.format(volCheck)
#         print 'Adjusting weights by skew'
#         temp = []
#         for i,skew in enumerate(skews):
#             temp.append(IBZ.weights[i]*(skew)**(1/4.0))    
#         IBZ.weights = temp/self.ravg**3   
#         IBZ.weights = IBZ.weights/self.ravg**3 #to scale them to order(1).  


    def meshInitCubic(self,type,eps):
        '''Add a cubic mesh to the interior, . If any 2 or 3 of the facet planes are 
        orthogonal, align the cubic mesh with their normals.       
        Remove any points within self.dw*self.wallClose from any wall        '''
        a = 1.0
        cubicLVs = identity(3)
        #test facet points for orthogonality
        rs = []
        pairs = []
        triples = []
        uvecs = self.IBZ.bounds[0]
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
        if self.initSrch is not None:
            self.IBZ = self.searchNmax(cubicLVs,aKcubConv,sites) 
        else:
            shift = array([1,1,1])/8.0 * aKcubConv
            self.IBZ,nInside = self.fillMesh(cubicLVs,self.IBZ,shift,aKcubConv,sites)
        
        #Search over shift and rotation to find the most possible points inside

    def searchNmax(self,cubicLVs,aKcubConv,sites):
        '''Test an entire grid of init values'''
        cubicLVs0 = cubicLVs
        nShift = 5
        nTh = 10
        nPh = 20
        shiftDiv = 0.5*sqrt(3)/float(nShift)
        thDiv = 90/float(nTh) #deg
        phDiv = 180/float(nPh)
        shifts = [i*shiftDiv*self.ravg*array([1,1,1]) for i in range(nShift)]
        thetas = [i*thDiv for i in range(nTh)]
        phis = [i*phDiv for i in range(nPh)]
        nMax = 0
        closestLog = 10
        isearch = 0
        bestN = 0
        IBZ = deepcopy(self.IBZ)
        bestIBZ = deepcopy(self.IBZ)
        
        for shift in shifts:
#             print 'shift',shift
            for theta in thetas:
                for phi in phis:
                    isearch += 1
                    Rmat = dot(
                        array([[1,0,0],
                        [0,cos(theta),-sin(theta)],
                        [0, sin(theta), cos(theta)]]),
                        array([[0,0,1],
                        [cos(phi),-sin(phi),0],
                        [sin(phi), cos(phi),0]]) )
                        
                    cubicLVs = dot(Rmat,cubicLVs0)
                    for i, site in enumerate(deepcopy(sites)):
                        sites[i] = dot(Rmat,site)        
                    IBZ,nInside = self.fillMesh(cubicLVs,IBZ,shift,aKcubConv,sites)
#                     self.facetsMeshMathFile(IBZ,'IBZinit_{}'.format(isearch),None)
#                     print isearch,'theta,phi',theta,phi,'n',nInside
#                     print 'nInside', nInside
                    if nInside > nMax: 
                        nMax = nInside
                        if self.initSrch == 'max':
                            bestN = nInside
                            bestIBZ = deepcopy(IBZ)
                            besti = isearch
                        print 'Step {}: Nmax'.format(isearch),nMax, shift, theta, phi
                    if self.initSrch != 'max':
                        closeLog = abs(log(nInside/(self.searchInitFactor*self.nTargetIBZ)))
                        if closeLog < closestLog:
                            closestLog = closeLog
                            bestN = nInside
                            bestIBZ = deepcopy(IBZ)
                        
#                         print 'Found nInside {} vs target N {} after {} steps'\
#                         .format(nInside,self.searchInitFactor*self.nTargetIBZ,isearch)
#                         return IBZ
        
#         else:
# #             print 'Did not find nInside within {}% of target N {}'.format(tol*100,self.nTargetIBZ)
#             print 'nInside {} vs target N {} after searching entire grid'.format(bestN,self.nTargetIBZ)
# #             print 'nInside {} is less than target N {}'.format(bestN,self.nTargetIBZ)
#              return bestIBZ
        if self.initSrch == 'max':
            print 'Maximum nInside {}(step {}) found vs target N {}'.format(bestN,besti,self.nTargetIBZ)
        else:
            print 'nInside {} is closest to adjusted target N {}'.format(bestN,self.searchInitFactor*self.nTargetIBZ)
        return bestIBZ
        
    def fillMesh(self,cubicLVs,IBZ,shift,aKcubConv,sites):
        #Find the extremes in each cubLV direction:
        intMaxs = [] #factors of aKcubConv
        intMins = []
        for i in range(3):
            projs = []
            for point in IBZ.fpoints:
                shifted = point + shift
                projs.append(dot(cubicLVs[:,i],shifted)/aKcubConv**2)
            intMaxs.append(int(ceil(max(projs)))+2) #optimize: Why is +2 required with shift of 1/2,1/2,1/2 on cubic?
            intMins.append(int(floor(min(projs)))-1)#optimize: Is -1 required?       
        #Create the cubic mesh inside the irreducible BZ
        IBZ.mesh = []
        IBZ.weights = []
        weightsInside = 0
        nInside = 0         
        ik = 0       
        for i in range(intMins[0],intMaxs[0]):
            for j in range(intMins[1],intMaxs[1]):
                for k in range(intMins[2],intMaxs[2]):
                    lvec = i*cubicLVs[:,0]+j*cubicLVs[:,1]+k*cubicLVs[:,2]
                    for iS, site in enumerate(sites):
                        ik+=1
                        kpoint = lvec + shift + site
                        if isInside(kpoint,IBZ.bounds,self.dw*self.wallClose):  #Can't be closer than self.dw*self.wallClose to a wall
                            nInside += 1
                            IBZ.mesh.append(kpoint) 
#         print 'Points inside', nInside
        return IBZ,nInside

    def dynamic(self,eps):
        ''' '''
        self.facetsMeshMathFile(self.IBZ,'IBZmeshInit',None)
        self.relax()
        self.facetsMeshMathFile(self.IBZ,'IBZmesh',None)
        return

    def relax(self):
        '''Minimization of the potential energy.
        The energy must be a function of 1-D inputs, so we flatten the points into components '''
        
        epsilon = self.ravg/100
        comps = array(self.IBZ.mesh).flatten()
        energy = self.minSteepest(comps,self.eps) 
        return

    def minSteepest(self,x0,eps):
        ''' Steepest descent works better in tests than sophisticated methods such as 
        conjugate gradient or using a Hessian method (see intro to P. Pulay, Chem. Phys. Lett. 73, 393 (1980)). 
        We move in the direction of -gradient, starting with a default step, and decreasing it by 2 until
        a lower energy and lower norm(gradient) is found.   The default step is increased if the previous step succeeded without adjustment.
        If the step must be lowered to find a lower energy and lower norm(gradient), the default is lowered. This worked
        better than the Barzilai-Borwein method that tries to predict an optimal step.    
        
'''

        itermax = 100
        gnormTol = 0.001
        minstep = 0.000001
        xold = x0
        xnew = x0
        fold,gold = self.enerGrad(xold)
        gnew = gold
        fnew = fold
        gnormold = norm(gold)
        gnormnew = gnormold
        fstart = fold; gstart = gold; gnormstart = gnormold
        method = 'steepest'
#         xolder =  xold + 0.01*gold #go backwards
#         golder = dot(H,(xolder-xold))
#         folder = dot(gold,(xolder-xold))
        print 'energy_0',fold, 'gnorm',gnormold #,'grad', gold
        iIter = 0
        step = 1.0 #* minstep
        atMinStep = False
        while iIter < itermax and gnormold > gnormTol and not atMinStep:
            print iIter, #progress bar
            method = 'steepest'
            lower = False
            while not lower:
                if step < minstep:
                    print 'minimum step reached: {}'.format(step) 
                    atMinStep = True
                    break
                if method == 'steepest':
                    xnew = xold - step*gold
                currPoints = xnew.reshape((len(self.IBZ.mesh),3))
                inside = True
                for point in currPoints:
                    if not isInside(point,self.IBZ.bounds,eps):
                        print 'point is outside IBZ...reduce step size:',point
                        inside = False
                        break                   
                if inside:
                    fnew,gnew = self.enerGrad(xnew)
                    gnormnew = norm(gnew)
                    if fnew<fold and gnormnew < gnormold:
                        lower = True
                        self.IBZ.mesh = currPoints.tolist()
                step /= 2
            step *= 4
            xold = xnew
            fold = fnew
            gold = gnew
            gnormold = gnormnew
            iIter += 1                   
        self.IBZ.mesh = currPoints.tolist()
        print 'For {} points in IBZ and {} steps'.format(len(self.IBZ.mesh),iIter)
        print '\tStarting energy',fstart, 'gnorm',gnormstart
        print '\tEnding energy',fnew,'gnorm',gnormnew, 'step',step#, 'grad', gnew
        print 'Energy/N',fnew/len(self.IBZ.mesh)
        if gnormnew <= gnormTol:
            print '\nSuccess after {} iterations'.format(iIter)
        elif iIter == itermax:
            print '\nExceeded maximum number of iterations ({}), while gnorm {} is greater than the tolerance {}'.format(itermax,gnormnew,gnormTol)
        if not (fnew < fstart and gnormnew < gnormstart):
            sys.exit('Did not find a lower energy and force norm: stop')
        return fnew

    def enerGrad(self,comps):
        '''Returns the total energy, gradient (-forces), 
        using comps, which are the flattened components of the current positions  '''
#         print 'oldindvecs',self.oldindVecs
        self.forces = zeros((len(self.IBZ.mesh),3))
        self.wallForce = zeros(len(self.IBZ.facets))
#         self.wallPress = zeros(len(self.IBZ.facets))
        IBZvecs = comps.reshape((len(self.IBZ.mesh),3))
        p = self.power
        wp = self.wallPower
        wallfact = self.wallfactor
        interfact = self.interfactor
        etot = 0
        for i,ri in enumerate(IBZvecs):
            #wall forces
            for iw, u in enumerate(self.IBZ.bounds[0]):
                ro = self.IBZ.bounds[1][iw]
                d = ro-dot(ri,u)+ self.wallOffset*self.dw #distance from plane to ri offset factor allows points to move closer to walls. 
                if d<0:
                    print '\nri,ro,u, dot(ri,u),d'
                    print ri,ro,u, dot(ri,u), d 
                    sys.exit('Error. Point {} in enerGrad is not in the IBZ.'.format(i+1))
                fmag = wallfact*(d/self.dw)**(-wp)  #dimensionless
                etot += wallfact*self.dw/abs(-wp+1)*(d/self.dw)**(-wp+1)#Units of length. Both F and E can't be dimensionless unless we make the positions dimensionless.
                self.forces[i] += -u*fmag
                self.wallForce[iw] += fmag #since forces are normal to plane, we sum the magnitudesrce',-u*fmag,fmag
            #vertext pull forces
#             for vert in self.IBZ.fpoints:
#                 d = norm(ri-vert)
#                 self.forces[i] += -self.vertexPull*(d/self.df)**(-p)*(ri-vert)/d #pull not push
#                 etot +=  self.vertexPull*self.df/abs(-p+1)*(d/self.df)**(-p+1)
#            inter-point forces
            for j, rj in enumerate(IBZvecs):
                if i!=j:
                    d = norm(ri-rj)
#                     print 'Inter d,f', d,interfact*(d/self.df)**(-p)*(ri-rj)/d
                    self.forces[i] += interfact*(d/self.df)**(-p)*(ri-rj)/d
                    if j>i: #don't overcount
                        etot += interfact*self.df/abs(-p+1)*(d/self.df)**(-p+1)
#         for i,fac in enumerate(self.facets):
#             area = convexH(planar3dTo2d(fac,self.eps)).volume  # for 2d problems, the "volume" returned is the area, and the "area" is the perimeter
#             self.wallPress[i] = self.wallForce[i]/area
#         print 'Pressure avg', mean(self.wallPress)
        return etot, -self.forces.flatten() #gradient is opposite the forces.
        
    def writeKpoints(self):
        nk = len(self.IBZ.mesh)
        totw = sum(self.IBZ.weights)
        lines = []
        lines.append('Packing of IBZ (Bret Hess, BYU). Total weights: {:12.8f} \n'.format(totw))#(vs 1.0 per general point without symmetry\n'.format(totw))
        lines.append('{}\n'.format(nk))
        lines.append('Reciprocal\n') #direct coordinates!...vasp doesn't read cartesian kpoints right
        etot = 0
        for ik,kpoint in enumerate(self.IBZ.mesh):
            kpointDir = directFromCart(self.B,kpoint)
            lines.append('{:15.12f}  {:15.12f}  {:15.12f}  {:15.12f}\n'\
                         .format(kpointDir[0],kpointDir[1],kpointDir[2],self.IBZ.weights[ik]))
        writefile(lines,'KPOINTS') 
                
       
    def getNeighbors(self,kpoint,IBZ,eps):
        '''Search a sphere around the kpoint and collect neighbors.
        Then if outside the IBZ, move point into IBZ by symmetry.  Search another sphere.
        Return the neighbors        '''
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
        As a convention, choose to keep points and volumes that are closest to the 111 direction. 
        
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
 
        self.IBZvolCut = 1.0 
        oldBZ = deepcopy(BZ)
        oldIBZvolCut = copy(self.IBZvolCut)
        reflOps = []
        rotOps = []
        rotReflOps = []
#         inversion = False
        inversion = True
        for iop in range(self.nops):
            op = self.symops[:,:,iop]
            if areEqual(trace(op),-3.0,eps):
                inversionNatural = True  
#                 inversion = True  
#                 
            if areEqual(abs(trace(op)),3.0,eps):#skip identity
                continue
            evals,evecs = eig(op)
            evecs = array([evec for evec in evecs])
            if areEqual(det(op),-1.0,eps) and not allclose(imag(evecs),zeros((3,3)),atol=eps): #improper rotation
                rotOps.append(-op)
            elif areEqual(det(op),1.0,eps)  : #rotation
                rotOps.append(op)
            else: #includes
                reflOps.append(op)
        for op1 in rotOps: 
            for op2 in reflOps:
                rotReflOps.append(dot(op1,op2))
        #add inversion if needed:
#         if not inversion:
#             self.nops+=1
#             print 'Symmetry operations (inversion added):', self.nops
#             inversion = True
#         else:
#             print 'Symmetry operations:', self.nops 
        print 'Symmetry operations (no inversion added):', self.nops
        print '\nReflection ops:'
        for iop, op in enumerate(reflOps):    
            print '\nsymop',iop; print op; print                
            evals,evecs = eig(op)
            evecs = array([evec for evec in evecs])                          
            if len(where(areEqual(evals,-1.0,eps))[0] )> 1: 
                evals = -evals 
                print 'Improper reflection'
            evec = evecs[:,where(areEqual(evals,-1.0,eps))[0][0]]
            u1 = self.choose111(evec,eps) 
            print 'reflection u1',u1
            BZ = self.cutCell(u1,0.0,BZ,eps)
            self.testCuts(BZ,oldBZ,oldIBZvolCut,'refl_{}'.format(iop),eps)
            if areEqual(self.IBZvolCut,self.nops,1e-2): return BZ
        print '\nRotation ops:'
        for iop, op in enumerate(rotOps):
            print '\nsymop',iop;print op ;print                
            evals,evecs = eig(op)
            evecs = array([evec for evec in evecs])                          
            ievec = where(areEqual(evals,1.0,eps))[0][0]
            evec = evecs[:,ievec] #axis
            ds = []
            allPoints = BZ.fpoints
            for vec in allPoints:
                if areEqual(abs(dot(evec,vec)),norm(vec),eps): #axis and vec are parallel...don't want this one.
                     ds.append(100)
                else:
                    ds.append(norm(vec - evec*dot(evec,vec)))
            allPoints = [point for (d,point) in sorted(zip(ds,allPoints),key = lambda x: x[0])]#sort by distance
            #pnto = allPoints[0] #old method
            for ip,pnto in enumerate(allPoints):#
                print ip,
                pntp = dot(op,pnto)
                if not among(pntp,allPoints,eps):
                    break
                #the plane to cut is the plane of O and axis, so take normal perpendicular to vector O.                   )
                tempvec0 = cross(evec,pnto)
                if not allclose(tempvec0,0.0,eps):
                    u1 = self.choose111(tempvec0/norm(tempvec0),eps)                             
                         #the plane to cut is the plane of O and axis, so take normal perpendicular to vector O.                   )
                    tempvec = cross(evec,pnto)
                    u1 = self.choose111(tempvec/norm(tempvec),eps)
#                         print 'u1',u1
                    BZ = self.cutCell(u1,0.0,BZ,eps)
                    tempvec = cross(evec,pntp)/norm(cross(evec,pntp))#2nd cut plane for roation
                    if not allclose(tempvec,-u1,atol=eps): #don't cut again if this is a Pi rotation
                        if abs(dot(tempvec, array([1,1,1])))>eps:
                            u2 = self.choose111(tempvec/norm(tempvec),eps)
                        else:
                            u2 = -dot(op,u1)
                        BZ = self.cutCell(u2,0.0,BZ,eps)                        
                    if self.testCuts(BZ,oldBZ,oldIBZvolCut,'rot_{}'.format(iop),eps):
                        break 
            if areEqual(self.IBZvolCut,self.nops,1e-2): return BZ 
        if inversion and not areEqual(self.IBZvolCut,self.nops,eps):
            if not inversionNatural:
                print 'Inversion added'
                self.nops = 2*self.nops
            anyDups,point1,point2 = makesDups(array([[-1,0,0],[0,-1,0],[0,0,-1]]),BZ,eps)
            if anyDups:
                u1 = self.choose111(point1/norm(point1),eps)
                BZ = self.cutCell(u1,0.0,BZ,eps)
                print 'Inversion', u1
                self.testCuts(BZ,oldBZ,oldIBZvolCut,'inv',eps)
                BZ.fpoints = flatVecsList(cell.facets,eps)
            BZ.fpoints = flatVecsList(cell.facets,eps)
            return BZ           
        if not areEqual(self.IBZvolCut,self.nops,eps):
            print ('Fail: Volume not reduced by factor equal to the number of symmetry operations')
            BZ.fpoints = flatVecsList(cell.facets,eps)
            return BZ
        else:
            BZ.fpoints = flatVecsList(cell.facets,eps)
            return BZ
   
    def testCuts(self,BZ,oldBZ,oldIBZvolCut,tag,eps):
        getBoundsFacets(BZ,eps)
        BZ.fpoints = flatVecsList(BZ.facets,eps) 
        self.facetsMathFile(BZ,'{}'.format(tag))
        opDone = False               
        try:
            BZ.volume = convexH(BZ.fpoints).volume
            self.IBZvolCut = det(self.B)/BZ.volume
            print 'Cut vol BZ / Vol IBZ', self.IBZvolCut  
            if self.IBZvolCut != oldIBZvolCut and isinteger(self.IBZvolCut,1e-2) and self.IBZvolCut<=self.nops+eps:
                print 'OK'
                iopDone = True
                oldBZ = deepcopy(BZ)
                oldIBZvolCut = copy(self.IBZvolCut)
                opDone = True
            elif self.IBZvolCut>self.nops+eps:
                opDone = True
            else:
                print 'Noninteger or no change.'
                BZ = deepcopy(oldBZ)
                self.IBZvolCut = copy(oldIBZvolCut )        
        except:
            print 'Zero volume obtained.'
            BZ = deepcopy(oldBZ)
            self.IBZvolCut = copy(oldIBZvolCut)
        return opDone
    
    def facetsMathToStr(self,strOut,cell,label,axes = False, color = 'Red'):
        ''' Mathematica output for drawing the facets of a cell'''
        strOut += '{}'.format(label)+' = Graphics3D[{'+'{}'.format(color)+', Thick,{'
        for ifac, facet in enumerate(cell.facets):
            facet = list(trimSmall(array(facet)))
            strOut += 'Line[{'
            for point in facet:
                strOut += '{'+'{:12.8f},{:12.8f},{:12.8f}'.format(point[0],point[1],point[2])+'},'
            strOut += '{'+'{:12.8f},{:12.8f},{:12.8f}'.format(facet[0][0],facet[0][1],facet[0][2])+'}' #back to first
            strOut += '}]'
            if ifac < len(cell.facets)-1:
                strOut += ','
        strOut += '}}'
        if axes:
            strOut += ',Axes -> True,AxesLabel -> {"x", "y", "z"}]'
        else:
            strOut += ']'  
        return strOut
    
    def facetsMathFile(self,cell,tag):
        '''Output for Mathematica graphics drawing BZ facets'''
        strOut = ''
        strOut = self.facetsMathToStr(strOut,cell,'s','True','Red'); 
        strOut += ';\nShow[s]' 
#         strOut += '}];'
        writefile(strOut,'cell_{}.m'.format(tag))     
            
    def facetsMeshMathFile(self,cell,tag,color):
        '''Output for Mathematica graphics drawing BZ facets and spheres at each mesh point'''
        strOut = ''
        strOut = self.facetsMathToStr(strOut,cell,'s','True','Red'); 
        strOut += ';\n p=Graphics3D[{'
        cell.mesh = list(trimSmall(array(cell.mesh)))
        for ipoint,point in enumerate(cell.mesh):
            if color != None:
                strOut += '{},'.format(color)
            strOut += 'Opacity[0.7],Sphere[{' + '{:12.8f},{:12.8f},{:12.8f}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(self.rpacking)
            if ipoint < len(cell.mesh) -1:
                strOut += ','
        strOut += '}];\nShow[s,p]'
        writefile(strOut,'cell_{}.m'.format(tag))         
    def facetsMeshVCMathFile(self,BZ,meshFacets):
        '''Output for Mathematica graphics drawing the facets of each mesh point
        cell, as well as the borders of the BZ'''
        strOut = ''
        strOut = self.facetsMathToStr(strOut,BZ,'s','True','Red');  
        showCommand = 'Show[s,'  
        strOut += ';\n'
        for ipoint,point in enumerate(BZ.mesh):
            tCell = cell()
            tCell.facets = [[]]*len(meshFacets[ipoint])
            for ifac,facet in enumerate(meshFacets[ipoint]):
                facet = list(trimSmall(array(facet)))
                temp = []
                for facetPoint in facet:
                    temp.append(facetPoint + point)
                tCell.facets[ifac] = temp
            strOut = self.facetsMathToStr(strOut,tCell,'v{}'.format(ipoint),False,'RandomColor[]')
            strOut+=';'
            showCommand += 'v{}'.format(ipoint)
            if ipoint < len(BZ.mesh)-1:
                showCommand += ','
        showCommand += ']'
#         strOut+=';'
        strOut+=showCommand
        writefile(strOut,'IBZmeshFacets.m')

    def mathPrintPlanes(self,bounds):
        print 'r=RegionPlot3D[',
        for iu, uvec in enumerate(bounds[0]):
            print '{}x+{}y+{}z<={}'.format(uvec[0],uvec[1],uvec[2],bounds[1][iu]), #write plane equations for mathematica
            if iu < len(bounds[0])-1:
                print '&&'
        print ', {x, -2, 2}, {y, -2, 2}, {z, -2, 2}, PlotStyle -> Opacity[0.3]]\n'