from numpy import array, arccos, dot, pi,  floor, sum, sqrt, exp, log, matrix
from numpy.matlib import zeros #creates np.matrix rather than array
from numpy.linalg import norm, det
from numpy import int as np_int
from random import random, randrange
import os, subprocess, sys, time

'''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer matrix
So K = B Inv(M) 
 M   =   [a g h]
         [b c j]
         [d e f]  
         
Minimization scheme'''

def orthdef(latt): #columns as vectors
    od = norm(latt[:,0])*norm(latt[:,1])*norm(latt[:,2])/det(latt)
    return od

def changewhich(M,B):
    bestgrad = 0
    bestdel = zeros((3,3),dtype=np_int)
    Mold = M
    oldcost = cost(Mold,B)
#    print 'old cost', oldcost
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
        od = 10
    else:
        od = orthdef(K.vecs)
    return(od)

class lattice(object): #reciprocal lattice
    def __init__(self):         
        self.vecs = []
        self.det = []
        self.Nmesh = []
#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
Nmesh=640
M = zeros((3,3),dtype=np_int)
print 'Target N kpoints', Nmesh
a = int(Nmesh**(1/3.0)); c = a; f = int(Nmesh/a/c) 
print a,c,f #s
M[0,0]= a
M[1,1]= c
M[2,2]= f
alat = 2*sqrt(2)
clat = alat*4/3
B =  lattice()
B.vecs = matrix((
  [   -alat/2,  alat/2,   alat/2],
  [   alat/2,  -alat/2,   alat/2],
  [   clat/2,   clat/2,   -clat/2]
  ), dtype=float)
B.det = det(B.vecs)
B.Nmesh = Nmesh
print B.vecs
print 'Det', B.det
print 'Cost (Orth Defect)', orthdef(B.vecs)

#sys.exit('stop')
maxsteps = 1000
min_od = 100
istep = 0
while istep<maxsteps:
    print '\nstep %i  ' % istep
#    print 'Current K mesh\n', B.vecs*M.I
    bestindex = changewhich(M,B)
#    print bestindex
    if bestindex[2]==0:#found minimum
        break
    else:
        M[bestindex[0],bestindex[1]] += bestindex[2]
    newcost = cost(M,B)
#    print M
    print 'New cost: %f6.3\n' % newcost
    istep += 1
    
if istep < maxsteps:
    print 'Found minimum cost %f6.3 after %i steps' % (newcost,istep)
    print 'Best M'; print M
    K = lattice();K.vecs = B.vecs*M.I; K.det = det(K.vecs)
    print 'Best K mesh\n', K.vecs  
    print 'Number of mesh points',B.det/K.det
else:
    print 'Ended without minimum after maximum %i steps' % istep

            

            
#            for i in range(ntrials):
#                b = randrange(c)
#                d = randrange(f)
#                e = randrange(f)
#                print '    ', b,d,e

                        
            
            

#if __name__ == "__main__":
#        import sys
#        for index in xrange(1,len(sys.argv)):
#                print "Factors for %s : %s" %(sys.argv[index], str(factor(int(sys.argv[index]))))

#            for b in range(c):
#                hnf[1,0] = b
#                for d in range(f):
#                    hnf[2,0] = d
#                    for e in range(f):
#                         hnf[2,1] = e
#                         print hnf