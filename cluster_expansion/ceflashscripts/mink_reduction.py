#!/usr/bin/env python
from ctypes import byref, cdll, c_double
from numpy import dot, cross, abs, zeros, array
from numpy.linalg import norm
import sys
sys.path.append("../classes")
from poscar import POSCAR as P

rlines  = [i.strip() for i in sys.stdin.readlines()] #assumes use of  < POSCAR

pos = P(lines=rlines)

print "real  od:",pos.orthogonality_defect
print "recip od:",pos.rod
print "original lattice:"
print array(pos.avecs)
print
print "k-space lattice:"
print array(pos.bvecs)
print
print "Test by reducing real lattice:"
print pos.mink_reduce(1e-10)

print "Test by reducing k-space:"
print pos.kmink_reduce(1e-1)
# Note: the reducing algorithms do NOT change ANY values.  Have to 

pos.avecs = pos.mink_reduce(1e-10)
print 'New real lattice orth defect'  
print 'New recip lattice orth defect'  
print 'New avecs'
print pos.rod
print pos.orthogonality_defect
print pos.avecs