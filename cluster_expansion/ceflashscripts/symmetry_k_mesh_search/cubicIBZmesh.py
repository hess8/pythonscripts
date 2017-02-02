from meshpy.tet import MeshInfo, build, Options

import os, subprocess,sys,re,time
from numpy import (mod,dot,cross,transpose, rint,floor,ceil,zeros,array,sqrt,
                   average,std,amax,amin,int32,sort,count_nonzero,
                   delete,mean,square,argmax,argmin,insert,s_,concatenate)
# from scipy.optimize import fmin_cg
# from scipy.spatial import Delaunay as delaunay
from numpy.linalg import inv, norm, det
from numpy.random import rand
from copy import copy,deepcopy
from sched import scheduler
from itertools import chain
from matplotlib.pyplot import (subplots,savefig,imshow,close,plot,title,xlabel,
                               ylabel,figure,show,scatter,triplot,hist)
import matplotlib.image as mpimg
import datetime
from _ast import operator
sys.path.append('/bluehome2/bch/pythonscripts/hesslib/')
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/')
sys.path.append('/fslhome/bch/graphener/graphener')
# os.listdir('/home/hessb/research/pythonscriptsRep/pythonscripts/hesslib/')
# sys.path.append('/home/hessb/research/pythonscriptsRep/pythonscripts/cluster_expansion/ceflashscripts/')


# from conjGradMin.optimize import fmin_cg

# from kmeshroutines import (svmesh,svmeshNoCheck,svmesh1freedir, lattice_vecs, lattice, surfvol,
#     orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,
#     load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double,
#     getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection,intoVoronoi,intoCell,
#     reverseStructured,isInVoronoi)

def timestamp():
    return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())

def plotPoints(pts,tag,tets,mainDir):
    plot2dPts(pts,tag,'x','y',tets,mainDir) 
    plot2dPts(pts,tag,'x','z',tets,mainDir)

def plot2dPts(pts,tag,ax0,ax1,tets,mainDir):
    fig = figure()
    i0 = ('x','y','z').index(ax0)
    i1 = ('x','y','z').index(ax1)
    scatter(pts[:][i0],pts[:][i1])
    ind0 = []
    ind1 = []
    if len(tets)>0: 
        pairs = []
        for itet, tet in enumerate(tets):
            pairs.append([pts[tet[0]],pts[tet[1]]])
            pairs.append([pts[tet[0]],pts[tet[2]]])
            pairs.append([pts[tet[0]],pts[tet[3]]])
            pairs.append([pts[tet[1]],pts[tet[2]]])
            pairs.append([pts[tet[1]],pts[tet[3]]])
            pairs.append([pts[tet[2]],pts[tet[3]]])
        for pair in pairs:
            plot([pair[0][i0],pair[1][i0]],[pair[0][i1],pair[1][i1]]) #draw line for tet edge
    xlabel('{}'.format(ax0))
    ylabel('{}'.format(ax1))
    name = '{}_{}_{}'.format(tag,ax0,ax1)
    title(name)
    fig.savefig(mainDir+'/'+ name+'.pdf');
    os.system('cp {} {}'.format(mainDir+'/'+ name+'.pdf',mainDir+ '/latest{}-{}'.format(ax0,ax1)+'.pdf' ))
#     show()
   
    
def plothist(vols):
    # the histogram of the data
#     figure()
    n, bins, patches = hist(vols, 50, facecolor='green')
    xlabel('Volume')
    ylabel('Number')
    show()
#     close()

def volTet2(vertPoints): 
        return abs(dot((vertPoints[0]-vertPoints[3]),cross((vertPoints[1]-vertPoints[3]),(vertPoints[2]-vertPoints[3])))/6.0)
 
# mainDir = '/home/hessb/Downloads/temp'
mainDir = '/fslhome/bch/trash'

os.chdir(mainDir)

mesh_info = MeshInfo()
# Cubic IBZ
# mesh_info.set_points([
#     (0,0,0), (0.5,0.5,0), (0,0.5,0), (0.5,0.5,0.5)
#     ])
# mesh_info.set_facets([
#     [0,1,2],
#     [0,1,3],
#     [1,2,3],
#     [3,2,0]
#     ])

#Cubic entire cubic voronoi cell....
a = 0.5; b = -0.5
mesh_info.set_points([
    (b,b,b), (b,b,a), (b,a,a),(b,a,b),(a,b,b), (a,b,a), (a,a,a),(a,a,b)
    ])
mesh_info.set_facets([
    [0,1,2,3],
    [0,1,5,4],
    [0,4,7,3],
    [6,2,1,5],
    [6,2,3,7],
    [6,7,4,5]])

opts = Options("VO9pa0.04") # Overriding 'pq'  

mesh = build(mesh_info,options=opts)
# mesh = build(mesh_info,max_volume = 0.1)
vols = []
tets = []
points = []
print "Mesh Points:"
for i, p in enumerate(mesh.points):
    points.append(array(mesh.points[i]))
    print i, p
# print "Point numbers in tetrahedra:"
for i, tet in enumerate(mesh.elements):
    tets.append(tet)
    vpoints = array([mesh.points[iv] for iv in tet])
    vols.append(volTet2(vpoints))
vols = array(vols)
vols.sort()
print 'Number of tets',len(vols)
print 'Average vol', mean(vols)
print 'Max vol', amax(vols)
print 'Min vol', amin(vols)
print 'Std dev', std(vols)
print 'Std/Avg', std(vols)/mean(vols)
# mesh.write_vtk("test.vtk")
#!/usr/bin/env python
# import numpy as np
# import matplotlib.mlab as mlab
# import matplotlib.pyplot as plt

plothist(vols)

plotPoints(points,timestamp(),tets,mainDir)

print