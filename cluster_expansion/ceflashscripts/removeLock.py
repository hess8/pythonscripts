#!/usr/bin/python
import time, os, subprocess

mainDir = '/fslhome/bch/cluster_expansion/hexagonal/aflow2x2adatoms/'
toRunFile = mainDir + 'jobs2run'
jobsfile = open(toRunFile,'r')
lines1 = jobsfile.readlines()
jobsfile.close()
for path in lines1:
    lockfile = path.replace('aflow.in','LOCK' )   
    os.system('rm '+lockfile)   
print 'done'