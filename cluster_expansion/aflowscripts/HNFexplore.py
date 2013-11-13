from numpy import array, arccos, dot, pi, zeros, floor, sum, sqrt
from numpy.linalg import norm
from numpy import int as np_int
from random import random, randrange
import os, subprocess, sys, time

def FindAllDivisors(x):
    divList = []
    y = 1
    while y <= sqrt(x):
        if x % y == 0:
            divList.append(y)
            divList.append(int(x / y))
        y += 1
    return sorted(divList)

natoms = 3
nkppra = 10000
nk = int(nkppra/natoms)
print 'N kpoints', nk
aspace = FindAllDivisors(nk)
cspace = aspace
nacf = len(aspace)*len(cspace) # of acf combinations
print 'Number of a,c,f combinations', nacf
time.sleep(5)
hcf = zeros((3,3),dtype=np_int)
ntrials = 100 #random b,d,e combinations tested, for each acf combination
itrial = 0
for a in aspace:
    hcf[0,0] = a
    for c in cspace:
        hcf[1,1] = c
        f = nk/a/c
        print a,c,f    
        if f>0:
            hcf[2,2] = f
            nposs = (c)*(f)*(f)# number of bde combinations
            
            frac = nposs/ntrials
            print frac
            

            
            for i in range(ntrials):
                b = randrange(c)
                d = randrange(f)
                e = randrange(f)
                print '    ', b,d,e

                        
            
            

#if __name__ == "__main__":
#        import sys
#        for index in xrange(1,len(sys.argv)):
#                print "Factors for %s : %s" %(sys.argv[index], str(factor(int(sys.argv[index]))))

#            for b in range(c):
#                hcf[1,0] = b
#                for d in range(f):
#                    hcf[2,0] = d
#                    for e in range(f):
#                         hcf[2,1] = e
#                         print hcf