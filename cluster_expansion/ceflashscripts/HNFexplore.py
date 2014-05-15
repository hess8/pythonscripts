from numpy import array, arccos, dot, pi,  floor, sum, sqrt, exp, log, matrix
from numpy.matlib import zeros
from numpy.linalg import norm, det
from numpy import int as np_int
from random import random, randrange
import os, subprocess, sys, time

'''The kmesh can be related to the reciprocal lattice B by  B = KH, where H is Hermite normal form (HNF) 
(or B = kM, where M is a different type of integer matrix). So K = B Inv(H) 
 HNF =   [a 0 0]
         [b c 0]
         [d e f]   Use positive integers, and b<c, d<f, e<f'''

def FindAllDivisors(x):
    divList = []
    y = 1
    while y <= sqrt(x):
        if x % y == 0:
            divList.append(y)
            divList.append(int(x / y))
        y += 1
    return sorted(divList)

def setRange(int1,ntrysmall, ntrylarge):
    '''if int1 is small, take all values less than that; 
    if it's large, take all values less than ntry, then distribute ntry more on a log basis, ''' 
    if int1 <= ntrysmall:
        tryrange = range(int1)
    else:
        tryrange = range(ntrysmall)
        for i in range(min([ntrylarge,int1-ntrysmall])):
            tryrange.append(ntrysmall+int(exp(random()*(log(int1)-log(ntrysmall)))))
    return tryrange

def orthdef(latt): #columns as vectors
    od = norm(latt[:,0])*norm(latt[:,1])*norm(latt[:,2])/det(latt)
    return od
#
#natoms = 3
#nkppra = 10000
#nk = int(nkppra/natoms)
nk=640
print 'N kpoints', nk
aspace = FindAllDivisors(nk)
print aspace
cspace = aspace
nacf = len(aspace)*len(cspace) # of acf combinations
print 'Number of a,c,f combinations', nacf
#time.sleep(5)
hnf = zeros((3,3),dtype=np_int)
ntrysmall = 12 #random b,d,e combinations tested, for each acf combination
ntrylarge = 0
itrial = 0
alat = 2*sqrt(2)
clat = alat*4/3
Blatt = matrix((
  [   -alat/2,  alat/2,   alat/2],
  [   alat/2,  -alat/2,   alat/2],
  [   clat/2,   clat/2,   -clat/2]
  ), dtype=float)
print Blatt
print 'Orth Defect', orthdef(Blatt)
#sys.exit('stop')
min_od = 100

for a in aspace:
    hnf[0,0] = a
    for c in cspace:
        hnf[1,1] = c
        f = nk/a/c
        if f>0:
            print 
            print  a,'a'
            print c,f,'c,f'   
            hnf[2,2] = f
            nposs = (c)*(f)*(f)# number of bde combinations
            brange = setRange(c,ntrysmall, ntrylarge)
            drange = setRange(f,ntrysmall, ntrylarge)
            erange = setRange(f,ntrysmall, ntrylarge)
            print max(brange), max(drange),max(erange)
            print brange
            print drange
            print erange
            for b in brange:
                hnf[1,0] = b
                for d in drange:
                    hnf[2,0] = d
                    for e in erange:
                        hnf[2,1] = e
                        print hnf
                        K = Blatt*hnf.I
#                        print 'Kmesh'; print K
#                        print 'det Kmesh',det(K)
                        od = orthdef(K)
                        print 'Orth Defect of mesh',  od
                        if od < min_od: min_od = od
                        time.sleep(2)
#                        sys.exit('stop')

            itrial = itrial + len(brange)*len(drange)*len(erange)
print'total trials', itrial
print 'Minimum orth defect', min_od
            

            
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