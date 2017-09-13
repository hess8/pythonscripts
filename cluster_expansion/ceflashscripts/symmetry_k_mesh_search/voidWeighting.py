'''
This version is passed parameters from poscar2mesh7, which is called by dynPackTest 
for searching method parameter space.
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

def addPlane(u,mag,planeList,eps):
    '''Add plane to a list that is of the form [[u's],[mags]]'''
    for ip,u2 in enumerate(planeList[0]):
        if allclose(u,u2,eps) and areEqual(mag,planeList[1][ip],eps):
            break
    else:
        planeList[0].append(u)
        planeList[1].append(mag)

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

def threePlaneIntersect(uRows,mags):  
    '''This routine is for three planes that will intersect at only one point, 
    as borders of facets, as we have for the Voronoi cell boundaries. 
    Planes are given by normals to ro, at the point ro.  All points r in the
    plane obey dot(r,ro) = dot(ro,ro) = ro^2 
      If they intersect, then 
        inv([[xo,yo,zo]
        [x1,y1,z1]
        [x2,y2,z2]] }  is not singular
        
    Or, more useful when we have some ro's that are zero:
      dot(r,uo) = ro defines points on the plane
      If they intersect, then 
        [[uxo,uyo,uzo]
        [ux1,uy1,uz1]
        [ux2,uy2,uz2]] } * [rx,ry,yz] = [ro1,ro2,ro3]
        
        Then the intersection point is
        
        [rx,ry,rz]  = 
      
        inv([[uxo,uyo,uzo]
        [ux1,uy1,uz1]
        [ux2,uy2,uz2]] } * (ro,r1,r2).  
        This should work even when ro = r1 = r2 =0, 
        which returns (0,0,0), the origin
      
    '''
    try:
        invuMat = inv(array(uRows))
        invOK = True
    except:
        invOK = False
        return None
    if invOK:
        point = trimSmall(dot(invuMat,mags))  
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
    allPlanes = [[cell.bounds[0][i] for i in range(len(cell.bounds[0]))], [cell.bounds[1][i] for i in range(len(cell.bounds[0]))]]
    allVerts = deepcopy(cell.fpoints)
    for ig  in grp:
        bndsLabels.append(ig)
    pairs = list(combinations(bndsLabels,2))
    for ig in grp:
        for pair in pairs:
            planesu3 = [boundVecs[ig]['uvec']]
            planesmag3 = [boundVecs[ig]['mag']]
            if not ig in pair:
                planesu3.append(boundVecs[pair[0]]['uvec'])
                planesmag3.append(boundVecs[pair[0]]['mag'])
                planesu3.append(boundVecs[pair[1]]['uvec'])
                planesmag3.append(boundVecs[pair[1]]['mag'])
                intersPt = threePlaneIntersect(planesu3,planesmag3) 
                if not intersPt is None:
                    if not isOutside(intersPt,cell.bounds,eps)\
                        and not among(intersPt,allVerts,eps):
                            allVerts = addVec(intersPt,allVerts,eps)
                            addPlane(planesu3[0],planesmag3[0],allPlanes,eps)
                            addPlane(planesu3[1],planesmag3[1],allPlanes,eps)
                            addPlane(planesu3[2],planesmag3[2],allPlanes,eps)                                 
    if len(allVerts)>0:
        #keep only the vertices that can be reached without crossing any plane
        newVerts = []
        for vert in allVerts:
            for ip,u in enumerate(allPlanes[0]):
                if dot(vert,u)>allPlanes[1][ip]+eps: 
                    break
            else: 
                newVerts.append(vert)
        #Start new to remove planes that don't host a vertex.
        tryPlanes = deepcopy(allPlanes)
        cell.bounds = [[],[]]
        for ip,u in enumerate(tryPlanes[0]):
            for vert in newVerts:
                if areEqual(dot(vert,u),tryPlanes[1][ip],eps):
                    addPlane(u,tryPlanes[1][ip],cell.bounds,eps)         
                    break
        cell.fpoints = newVerts
#     self.mathPrintPlanes(allPlanes)
    if len(cell.fpoints) >= 3 and type in ['BZ','MP']:
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

def isInside(vec,bounds,eps,extra = 0):
    '''Inside means on opposite side of the plane vs its normal vector.  Extra shifts the test plane
    outward, for an expanded inside volume'''
    inside = zeros(len(bounds[0]),dtype = bool)
    for iplane, uvec in enumerate(bounds[0]):
        if dot(vec,uvec) - extra < bounds[1][iplane] - eps: #point is inside this plane
            inside[iplane] = True
    return all(inside)

def isOutside(vec,boundaries,eps):
    for iplane, uvec in enumerate(boundaries[0]): 
        pvec = uvec*boundaries[1][iplane]           
        if dot(vec,uvec) > boundaries[1][iplane] + eps: #point is outside this plane
            return True
    return False

# def isJustOutside(vec,boundaries,dmax,eps):
#     for iplane, uvec in enumerate(boundaries[0]): 
#         pvec = uvec*boundaries[1][iplane]           
#         if  boundaries[1][iplane] + eps < dot(vec,uvec) < boundaries[1][iplane] + dmax + eps: #point is outside this plane
#             return True
#     return False

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
        if type in ['BZ','MP'] and not checkNext: break
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
    braggVecs = zeros(5*5*5-1,dtype = [('uvec', '3float'),('mag', 'float')])  
    #Exclude (0,0,0) in array dimensions (-1)
    ipoint = 0
    for i in range(-2,3):
        for j in range(-2,3):
            for k in range(-2,3):
                if not (i==j==k==0):
                    vec = trimSmall(0.5*(i*LVs[:,0] + j*LVs[:,1] + k*LVs[:,2]))
#                     vec = trimSmall(0.5*dot(LVs,array([i,j,k])))
                    mag = norm(vec)
                    braggVecs[ipoint]['uvec'] = vec/mag
                    braggVecs[ipoint]['mag'] = mag
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
        self.facets = [] #points arranged in facets
        self.fpoints = [] #points as a set (but a list)
        self.volume = None
        self.volumes = []
        self.center = None #body center of cell
        self.mesh = [] #centers of each voronoi cell, or kpoints
        self.weights = None
        self.vorVols = []
        self.details = None
    
class voidWeight(): 
    ''''''
    from numpy import zeros,array,mod
    from numpy.random import rand, uniform
        
    def pack(self,A,B,totatoms,aTypes,postype,aPos,targetNmesh,meshtype,path,params):
        paramLabels = ['']
        [symopsList, fracsList] = get_spaceGroup(transpose(A),aTypes,transpose(aPos),1e-3,postype.lower()[0] == 'd')
        self.nops = len(symopsList)
        self.symops = zeros((3,3,self.nops),dtype = float)
        for iop in range(len(symopsList)):
            self.symops[:,:,iop] = trimSmall(array(symopsList[iop]))
        self.B = B
#         print '\nB (Recip lattice vectors as columns',B
#         print 'method',method
        vol = abs(det(B))
        IBZvol = vol/float(self.nops)
#         self.ravg = (vol/targetNmesh)**(1/3.0) #distance if mesh were cubic. 
        self.ravg = (IBZvol/targetNmesh)**(1/3.0) #distance if mesh were cubic. 
        self.wallClose = float(params[0]) #
        self.initSrch = 'max'
        eps = self.ravg/300
        self.eps = eps
#         self.initSrch = None
#         self.initSrch = 'target'
        self.nTargetIBZ = targetNmesh; print 'targets are for IBZ, not full BZ'
        self.path = path      
        self.BZ = cell() #instance
        self.BZ.volume = vol
        braggVecs = getBraggVecs(self.B)
        self.BZ = getVorCell(braggVecs,self.BZ,'BZ',eps)
        self.facetsMathFile(self.BZ,'BZ') 
        self.IBZ = self.getIBZ(deepcopy(self.BZ),eps) #now irreducible BZ
        self.writeBounds()
#         self.facetsMathFile(self.IBZ,'IBZ') 
        self.meshInitCubic(meshtype,eps)
        self.facetsMeshMathFile(self.IBZ,self.IBZ.mesh,'IBZmesh',None)
        nKmax = 150
        print 'Limiting nK to {}'.format(nKmax)
        if 2 < len(self.IBZ.mesh) <= nKmax:
            OK = True
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
        '''
        
        xMake a standard voronoi cell for each mesh point. 
        xFind the max distance rmaxMP between a MP cell center an a vertex. 
        xIf a MP has a facet outside the IBZ bounds, it is cut.
        xFind the voids near the surface (min d_planes < rmaxMP).
        xUse sym and translation to find all IBZ points that have partners outside
        but are close to the surface. 
        void weights are distributed among IBZ points that are close to it.   
        
        Etot = sum(E_i*vol_i)/volIBZ  =  [sum(E_IBZMP_m*vol_m) + sum(E_void_q*vol_q)]/volIBZ.
        But we assume that each E_void can be interpolated from the calculated E_IBZMP's.  
        
        ****So we use some kind of interpolation...*** Linear at first.  
        
        
        E_void = sum(E_MPclose_i*
        '''
        allMPfacets = []
        surfPoints = []
        self.IBZ.weights = []
        for ip,point in enumerate(self.IBZ.mesh):
            print ip,
            ibzMP = self.prepMP(point)
            for fpoint in ibzMP.fpoints:
                cut = False
                for iplane, uvec in enumerate(self.IBZ.bounds[0]):
                    ro = self.IBZ.bounds[1][iplane]
                    d = ro - dot(uvec,fpoint)
                    if d < self.rmaxMP:
                        cut = True
                        surfPoints.append(point)                                      
                        ibzMP = self.cutCell(uvec,ro,ibzMP,eps) # we always keep the part that is "inside", opposite u
#             if cut:
#                 allMPfacets.append(ibzMP.facets)
            allMPfacets.append(ibzMP.facets)
            self.IBZ.vorVols.append(ibzMP.volume)               
#         self.facetsMeshVCMathFile(self.IBZ,allMPfacets,'IBZMesh')
        vMPs = sum(self.IBZ.vorVols)
        stdev = std(self.IBZ.vorVols)
        meanV = mean(self.IBZ.vorVols)
        volCheck = 0.01
        volDiff = vMPs - self.IBZ.volume        
        volDiffRel = volDiff/self.IBZ.volume
        print 'Total volume of point Vor cells',vMPs,'vs IBZ volume', self.IBZ.volume
        print 'Relative volume in MP Vor cells', volDiffRel,'Abs volume difference', volDiff, 'Std dev/mean',stdev/meanV
        self.IBZ.weights = self.IBZ.vorVols
        
        
        #find void centers, which are portions of points on the original packing lattice
#         that lie outside the IBZ 
        self.voids = cell()
        self.voids.volume = 0.0
        for io, point in enumerate(self.outPoints):
            ds = [norm(point-surfpoint) for surfpoint in surfPoints]
            if min(ds) < 2*self.rmaxMP: #candidate to host an IBZ void: it's "just outside"
                joMP = self.prepMP(point) # point just outside the IBZ on the mesh lattice
                for fpoint in joMP.fpoints:
                    for iplane, uvec in enumerate(self.IBZ.bounds[0]):
                        ro = self.IBZ.bounds[1][iplane]
                        d = ro - dot(uvec,fpoint)
                        if d < self.rmaxMP:                                    
                            joMP = self.cutCell(uvec,ro,joMP,eps) # we always keep the part that is "inside", opposite u                 
#                             print 'Surf point vol', point, joMP.volume
                if joMP.volume > self.eps**3:
                    self.voids.facets.append(joMP.facets)
                    self.voids.volume += joMP.volume
                    self.voids.volumes.append(joMP.volume)
                    self.voids.mesh.append(joMP.center)
        print 'Total volume Vor cells plus voids:',vMPs+self.voids.volume,'vs IBZ volume', self.IBZ.volume 
        self.facetsMeshVCMathFile(self.IBZ,allMPfacets,'IBZmesh') 
        self.facetsMeshVCMathFile(self.IBZ,self.voids.facets,'voids')                     
        if not areEqual(vMPs+self.voids.volume, self.IBZ.volume, volCheck*self.IBZ.volume):
            sys.exit('Stop: point Voronoi cells plus voids do not sum to the IBZ volume.') 
        #find distances of each void point to IBZ mesh points and their symmetry partners.    
        #find symmetry parterns (expandedMesh), which includes themselves
        rCutoff = 3.0*self.rpacking
        expandedMesh = [[]]*len(self.IBZ.mesh)
        for iIBZ,point in enumerate(self.IBZ.mesh):
            BZpartners = []
            for iop in range(self.nops):
                op = self.symops[:,:,iop]
                symPoint = dot(op,point)
                BZpartners.append(symPoint)  #includes IBZ points
            allPartners = deepcopy(BZpartners)
            for BZpoint in BZpartners:
                for i in [-1,0,1]:
                    for j in [-1,0,1]:
                        for k in [-1,0,1]:
                            if not i==j==k==0:
                                transPoint = BZpoint + i*self.B[0] + j*self.B[1] + k*self.B[2]
        #                         print 'transpoint',transPoint
                                if isInside(transPoint,self.IBZ.bounds,self.eps,rCutoff):
                                    allPartners.append(transPoint)
            expandedMesh[iIBZ]= deepcopy(allPartners)            
            self.facetsMeshMathFile(self.IBZ,expandedMesh[iIBZ],'expmesh_{}'.format(iIBZ),None)
        # Divide the volume of each void and add it to the volume of each mesh point, 
        # according to how close expandedMesh points (that are partners of the mesh point) 
        # is to the void point      
        for iv, vpoint in enumerate(self.voids.mesh):
            dweights = zeros(len(self.IBZ.mesh),dtype = float)
            for iIBZ in range(len(self.IBZ.mesh)):
                ds = [norm(epoint-vpoint) for epoint in expandedMesh[iIBZ]]
                print 'ds',sorted(ds)
                dmin = min(ds)
                imin = argmin(ds)
                if dmin > self.eps:
                    for d in ds:
                        dweights[iIBZ] += dmin/d # this keeps weights all at one or below
                else: #one point gets all the weight
                    dweights = zeros(len(self.IBZ.mesh),dtype = float)
                    dweights[iIBZ] = 1.0
                    break                
            dweights = dweights/sum(dweights) #now weights sum to 1
            #Divide volume in void:
            for iIBZ in range(len(self.IBZ.mesh)):
                self.IBZ.weights[iIBZ] += self.voids.volumes[iv] * dweights[iIBZ]             
        wtot = sum(self.IBZ.weights)
        print 'Total volume in reweighted IBZ MPs:',wtot,'vs IBZ volume', self.IBZ.volume                       
        if not areEqual(wtot, self.IBZ.volume, volCheck*self.IBZ.volume):
            sys.exit('Stop: point Voronoi cells plus voids do not sum to the IBZ volume.') 
        #normalize the seights so that a full interior point gets weight nops of symmetry
       
        
        self.IBZ.weights =  array(self.IBZ.weights)/self.MP.volume * self.nops
        print 'Weights:'
        for i, weight in enumerate(self.IBZ.weights):
            print i, weight
        print 'Sum', sum(self.IBZ.weights)                
        return

    def prepMP(self,kpoint):
        cutMP = deepcopy(self.MP)
        for ifac, facet in enumerate(self.MP.facets):
            temp = []
            for point in facet:
                temp.append(point + kpoint)
            cutMP.facets[ifac] = temp
        cutMP.fpoints = []
        for ipoint, point in enumerate(self.MP.fpoints):
            cutMP.fpoints.append(point + kpoint)
        cutMP.center = kpoint
        return cutMP

    def meshInitCubic(self,type,eps):
        '''Add a cubic mesh to the interior, . If any 2 or 3 of the facet planes are 
        orthogonal, align the cubic mesh with their normals.       
        Remove any points within self.dw*self.wallClose from any wall        '''
        a = 1.0
        cubicLVs = identity(3)
        if type == 'fcc':    
            volKcubConv = self.IBZ.volume/self.nTargetIBZ*4
            aKcubConv = volKcubConv**(1/3.0)
            cubicLVs = cubicLVs * aKcubConv
            sites = [array([0, 0 , 0]), 1/2.0*(cubicLVs[:,1]+cubicLVs[:,2]),\
                     1/2.0*(cubicLVs[:,0]+cubicLVs[:,2]), 1/2.0*(cubicLVs[:,0]+cubicLVs[:,1])]
            self.rpacking = 1/2.0/sqrt(2)*aKcubConv
            pf = 4*4/3.0*pi*(1/2.0/sqrt(2))**3  #0.74
        elif type == 'bcc':
            volKcubConv = self.IBZ.volume/self.nTargetIBZ*2
            aKcubConv = volKcubConv**(1/3.0)
            cubicLVs = cubicLVs * aKcubConv
            sites = [array([0, 0 , 0]), 1/2.0*(cubicLVs[:,0]+cubicLVs[:,1]+cubicLVs[:,2])]
            self.rpacking = sqrt(3)/4.0*aKcubConv
            pf = 2*4/3.0*pi*(sqrt(3)/4.0)**3 #0.68
        elif type == 'cub':
            volKcubConv = self.IBZ.volume/self.nTargetIBZ
            aKcubConv = volKcubConv**(1/3.0)
            cubicLVs = cubicLVs * aKcubConv
            sites = [array([0, 0 , 0])]
            self.rpacking = aKcubConv/2
            pf = 4/3.0*pi*(1/2.0)**3 #0.52
        else:
            sys.exit('stop. Type error in meshCubic.')

        if self.initSrch is not None:
            self.IBZ,self.outPoints = self.searchNmax(cubicLVs,aKcubConv,sites) 
            [shift,theta,phi] = self.IBZ.details
            Rmat = dot(
                array([[1,0,0], [0,cos(theta),-sin(theta)],[0, sin(theta), cos(theta)]]),
                array([[cos(phi),-sin(phi),0],[sin(phi), cos(phi),0],[0,0,1],]) )
            cubicLVs = dot(Rmat,cubicLVs)
            if type == 'fcc':    
                sites = [array([0, 0 , 0]), 1/2.0*(cubicLVs[:,1]+cubicLVs[:,2]),\
                         1/2.0*(cubicLVs[:,0]+cubicLVs[:,2]), 1/2.0*(cubicLVs[:,0]+cubicLVs[:,1])]
                self.meshPrimLVs = transpose(array(sites[1:]))
            elif type == 'bcc':
                sites = [array([0, 0 , 0]), 1/2.0*(cubicLVs[:,0]+cubicLVs[:,1]+cubicLVs[:,2])]
                self.meshPrimLVs = transpose(array([array(-1/2.0*(cubicLVs[:,0]+cubicLVs[:,1]+cubicLVs[:,2])),\
                            array(1/2.0*(cubicLVs[:,0]-cubicLVs[:,1]+cubicLVs[:,2])),\
                            array(1/2.0*(cubicLVs[:,0]+cubicLVs[:,1]-cubicLVs[:,2]))]))
            elif type == 'cub':
                sites = [array([0, 0 , 0])]
                self.meshPrimLVs = cubicLVs        
        else:
            shift = array([1,1,1])/8.0 * aKcubConv
            self.IBZ,self.outPoints = self.fillMesh(cubicLVs,self.IBZ,shift,aKcubConv,sites)


        MPbraggVecs = getBraggVecs(self.meshPrimLVs)
        self.MP = cell()
        self.MP.volume = self.IBZ.volume/self.nTargetIBZ
        self.MP = getVorCell(MPbraggVecs,self.MP,'MP',eps)
        self.rmaxMP = max([norm(point) for point in self.MP.fpoints])
        self.Vsphere = 4/3.0*pi*self.rpacking**3        
        #Search over shift and rotation to find the most possible points inside

    def searchNmax(self,cubicLVs,aKcubConv,sites):
        '''Test an entire grid of init values'''
        cubicLVs0 = cubicLVs
        nShift = 5
        
#         nTh = 10
#         nPh = 20
        
               
        nTh = 3
        nPh = 3
# 


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
        sites0 = deepcopy(sites)
        
        for shift in shifts:
#             print 'shift',shift
            for theta in thetas:
                for phi in phis:
                    isearch += 1
                    Rmat = dot(
                        array([[1,0,0], [0,cos(theta),-sin(theta)],[0, sin(theta), cos(theta)]]),
                        array([[cos(phi),-sin(phi),0],[sin(phi), cos(phi),0],[0,0,1],]) )
                      
                    cubicLVs = dot(Rmat,cubicLVs0)
                    for i, site in enumerate(sites0):
                        sites[i] = dot(Rmat,site)        
                    IBZ,nInside,outPoints = self.fillMesh(cubicLVs,IBZ,dot(Rmat,shift),aKcubConv,sites)
#                     self.facetsMeshMathFile(IBZ,IBZ.mesh'IBZinit_{}'.format(isearch),None)
#                     print isearch,'theta,phi',theta,phi,'n',nInside
#                     print 'nInside', nInside
                    if nInside > nMax: 
                        nMax = nInside
                        if self.initSrch == 'max':
                            bestN = nInside
                            bestOutside = outPoints
                            bestIBZ = deepcopy(IBZ)
                            bestIBZ.details = [shift,theta,phi]
                            besti = isearch
                        print 'Step {}: Nmax'.format(isearch),nMax, shift, theta, phi
                    if self.initSrch != 'max':
                        closeLog = abs(log(nInside/(self.nTargetIBZ)))
                        if closeLog < closestLog:
                            closestLog = closeLog
                            bestOutside = outPoints
                            bestN = nInside
                            bestIBZ = deepcopy(IBZ)
                            bestIBZ.details = [shift,theta,phi]
        if self.initSrch == 'max':
            print 'Maximum nInside {} (step {}) found vs target N {}'.format(bestN,besti,self.nTargetIBZ)
        else:
            print 'nInside {} is closest to adjusted target N {}'.format(bestN,self.nTargetIBZ)
        return bestIBZ,bestOutside
        
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
        outPoints = []
        nInside = 0         
        ik = 0       
        for i in range(intMins[0],intMaxs[0]):
            for j in range(intMins[1],intMaxs[1]):
                for k in range(intMins[2],intMaxs[2]):
                    lvec = i*cubicLVs[:,0]+j*cubicLVs[:,1]+k*cubicLVs[:,2]
                    for site in sites:
                        ik+=1
                        kpoint = lvec + shift + site
#                         if isOutside(kpoint,IBZ.bounds,self.eps):
#                             outPoints.append(kpoint)
#                         else:
#                             nInside += 1
#                             IBZ.mesh.append(kpoint) 
                        if isInside(kpoint,IBZ.bounds,self.eps,-self.wallClose*self.rpacking):
                            nInside += 1
                            IBZ.mesh.append(kpoint) 
                        else:
                            outPoints.append(kpoint)



                            
#         print 'Points inside', nInside
        return IBZ,nInside,outPoints

        
    def writeKpoints(self):
        nk = len(self.IBZ.mesh)
        totw = sum(self.IBZ.weights)
        lines = []
        lines.append('Packing of IBZ (Bret Hess, BYU). Total weights: {:12.8f} \n'.format(totw))#(vs 1.0 per general point without symmetry\n'.format(totw))
        lines.append('{}\n'.format(nk))
        lines.append('Reciprocal\n') #direct coordinates!...vasp doesn't read cartesian kpoints right
        for ik,kpoint in enumerate(self.IBZ.mesh):
            #direct coordinates!...vasp doesn't read cartesian kpoints right
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
                try:
                    cell.volume = convexH(cell.fpoints).volume    
                except:
                     cell.volume = 0                
        return cell      
#         else:
#             cell.volume = 0
#             return cell

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
            
    def facetsMeshMathFile(self,cell,points,tag,color):
        '''Output for Mathematica graphics drawing BZ facets and spheres at each mesh point'''
        strOut = ''
        strOut = self.facetsMathToStr(strOut,cell,'s','True','Red'); 
        strOut += ';\n p=Graphics3D[{'
        list(trimSmall(array(points)))
        for ipoint,point in enumerate(list(trimSmall(array(points)))):
            if color != None:
                strOut += '{},'.format(color)
            strOut += 'Opacity[0.7],Sphere[{' + '{:12.8f},{:12.8f},{:12.8f}'.format(point[0],point[1],point[2])+ '},'+'{}]'.format(self.rpacking)
            if ipoint < len(points) -1:
                strOut += ','
        strOut += '}];\nShow[s,p]'
        writefile(strOut,'cell_{}.m'.format(tag))   
              
    def facetsMeshVCMathFile(self,IBZ,allMeshFacets,tag):
        '''Output for Mathematica graphics drawing the facets of each mesh point
        cell, as well as the borders of the IBZ'''
        strOut = ''
        strOut = self.facetsMathToStr(strOut,IBZ,'s','True','Red');  
        showCommand = 'Show[s,'  
        strOut += ';\n'
        for i,facets in enumerate(allMeshFacets):
            tCell = cell()
            tCell.facets = facets
            strOut = self.facetsMathToStr(strOut,tCell,'v{}'.format(i),False,'RandomColor[]')
            strOut+=';'
            showCommand += 'v{}'.format(i)
            if i < len(allMeshFacets)-1:
                showCommand += ','
        showCommand += ']'
#         strOut+=';'
        strOut+=showCommand
        writefile(strOut,'facets_{}.m'.format(tag))

    def mathPrintPlanes(self,bounds):
        print 'r=RegionPlot3D[',
        for iu, uvec in enumerate(bounds[0]):
            print '{}x+{}y+{}z<={}'.format(uvec[0],uvec[1],uvec[2],bounds[1][iu]), #write plane equations for mathematica
            if iu < len(bounds[0])-1:
                print '&&'
        print ', {x, -2, 2}, {y, -2, 2}, {z, -2, 2}, PlotStyle -> Opacity[0.3]]\n'
            
    def writeBounds(self):
        lines = []
        for ib, u in enumerate(self.IBZ.bounds[0]):
            lines.append('{:12.8f} {:12.8f} {:12.8f} {:12.8f}\n'.format(u[0],u[1],u[2],self.IBZ.bounds[1][ib]))
        writefile(lines,'bounds')