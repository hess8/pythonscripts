import os, subprocess, sys, time

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice_vecs


from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, matrix, transpose,rint
from numpy.matlib import zeros #creates np.matrix rather than array
from numpy.linalg import norm, det
from numpy import int as np_int
from random import random, randrange


'''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer matrix
So K = B Inv(M) 
 M   =   [a g h]
         [b c j]
         [d e f]  
         
Minimization scheme'''

def orthdef(latt): #columns as vectors
    od = norm(latt[:,0])*norm(latt[:,1])*norm(latt[:,2])/det(latt)
    return od

def surfvol(vecs):
    vecs = transpose(vecs)
    u = norm(cross(vecs[0,:],vecs[1,:]))
    v = norm(cross(vecs[1,:],vecs[2,:]))
    w = norm(cross(vecs[2,:],vecs[0,:]))
    surfvol = 2*(u + v + w)/6/(det(vecs))**(2/3.0)
    return surfvol

def surfproduct(vecs):
    vecs = transpose(vecs)
    u = norm(cross(vecs[0,:],vecs[1,:]))
    v = norm(cross(vecs[1,:],vecs[2,:]))
    w = norm(cross(vecs[2,:],vecs[0,:]))
    surfproduct = sqrt(u * v * w)/det(vecs)
    return surfproduct
    
def changewhich(M,B):
    bestgrad = 0
    bestdel = zeros((3,3),dtype=np_int)
    Mold = M
    oldcost = cost(Mold,B)
    bestindex=[-1,-1,0]#initialize
    for i in range(3):
        for j in range(3):          
            M[i,j] += 1;delInc = cost(M,B)-oldcost; M[i,j] += -1           
            if delInc < 0 and delInc < bestgrad: bestindex = [i,j,1];bestgrad = delInc
            M[i,j] += -1;delDec = cost(M,B)-oldcost;M[i,j] += 1;
            if delDec < 0 and delDec < bestgrad: bestindex = [i,j,-1];bestgrad = delDec
#            print i,j,delInc,delDec
#            print
#            print i,j, best
    return bestindex

def cost(M,B):
    K = lattice()
    K.vecs = B.vecs*M.I;K.det = det(K.vecs)
    if B.det/K.det > 1.1*B.Nmesh or B.det/K.det < 0.9*B.Nmesh:
        cost = 10
    else:
#        cost = surfproduct(K.vecs)
        cost = surfvol(K.vecs)
#        cost = orthdef(K.vecs)
    return(cost)

class lattice(object): #reciprocal lattice
    def __init__(self):         
        self.vecs = []
        self.det = []
        self.Nmesh = []
#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=20000
M = zeros((3,3),dtype=np_int)
print 'Target N kpoints', Nmesh
a = rint(Nmesh**(1/3.0)); c = a; f = int(Nmesh/a/c) 
M[0,0]= a
M[1,1]= c
M[2,2]= f
print 'Starting M'
print M
B =  lattice()

##############BCT lattice
alat = 2*sqrt(2)
clat = alat*4/3
B.vecs = matrix((  
  [   -alat/2,  alat/2,   alat/2],
  [   alat/2,  -alat/2,   alat/2],
  [   clat/2,   clat/2,   -clat/2]
  ), dtype=float)
#############End BCT lattice

############## Any lattice
#angle1 = 90
#crystal = [1,1,1,angle1,angle1,angle1] #trigonal [a,b,c,alpha,beta,gamma]
#B.vecs = transpose(lattice_vecs(crystal))
#############End BCT lattice

B.det = det(B.vecs)
B.Nmesh = Nmesh
print B.vecs
print 'Det', B.det
print 'Orth Defect of B', orthdef(B.vecs)
print 'Surf/vol of B', surfvol(B.vecs)
print 'Surfproduct of B', surfproduct(B.vecs)
#sys.exit('stop')
maxsteps = 1000
min_od = 100
istep = 1
while istep<maxsteps:
#    print 'step %i  ' % istep
#    print 'Current K mesh\n', B.vecs*M.I
    bestindex = changewhich(M,B)
#    print bestindex
    if bestindex[2]==0:#found minimum
        newcost = cost(M,B)        
        break
    else:
        M[bestindex[0],bestindex[1]] += bestindex[2]
        newcost = cost(M,B)
    istep += 1
#    print M
#    print 'New cost: %f \n' % newcost
   
       
if istep < maxsteps:
    print
    print 'Found minimum after %i steps' % istep
    print 'Best M'; print M
    K = lattice();K.vecs = B.vecs*M.I; K.det = det(K.vecs)
    print 'Best K mesh\n', K.vecs  
    print 'Number of mesh points', B.det/K.det
    print 'Minimum cost', newcost
    print 'Orth defect',orthdef(K.vecs)
    print 'Surface/vol', surfvol(K.vecs)    
    print 'Mesh vector lengths', norm(K.vecs[:,0]),norm(K.vecs[:,1]),norm(K.vecs[:,2])

else:
    print 'Ended without minimum after maximum %i steps' % istep
