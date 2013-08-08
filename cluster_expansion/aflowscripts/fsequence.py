#!/usr/bin/python
import time, os, subprocess

mainDir = '/fslhome/bch/cluster_expansion/alir/'
toRunFile = mainDir + 'f1_50.dat'
file1 = open(toRunFile,'w')
for i in range(1,51):
    file1.write('%i  %i  f%s %s' % (i,i,i,'\n'))
file1.close() 
print 'done'