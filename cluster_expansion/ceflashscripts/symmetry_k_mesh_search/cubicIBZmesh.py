from meshpy.tet import MeshInfo, build

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


def volTet(vertPoints): 
        return dot((vertPoints[0]-vertPoints[3]),cross((vertPoints[1]-vertPoints[3]),(vertPoints[2]-vertPoints[3])))/6.0
 

mesh_info = MeshInfo()
mesh_info.set_points([
    (0,0,0), (0.25,0.25,0.25), (0.5,-0.5,0.5), (0,0,0.5)
    ])
mesh_info.set_facets([
    [0,1,2],
    [0,1,3],
    [1,2,3],
    [3,2,0]
    ])
mesh = build(mesh_info,max_volume = 0.01)
# mesh = build(mesh_info,max_volume = 0.1)
vols = []
print "Mesh Points:"
for i, p in enumerate(mesh.points):
    print i, p
print "Point numbers in tetrahedra:"
for i, tet in enumerate(mesh.elements):
    vpoints = [mesh.points[iv] for iv in tet]
    vols.append(volTet(vertPoints))
vols = array(vols)
print 'Average vol', mean(vols)
print 'Max vol', max(vols)
print 'Min vol', min(vols)
print 'Std dev', std(vols)
# mesh.write_vtk("test.vtk")