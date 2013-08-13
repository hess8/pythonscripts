#!/usr/bin/python
import time, os, subprocess

mainDir = '/fslhome/bch/cluster_expansion/alir/'
toRunFile = mainDir + 'f11000.dat'
file1 = open(toRunFile,'w')
for i in range(1,11001):
    file1.write('%i  %i  f%s %s' % (i,i,i,'\n'))
file1.close() 
print 'done'