import sys,os,subprocess,time,sys,shutil, numpy as np
#from ceScriptTools import runningJobs
dir = '/fslhome/bch/cluster_expansion/cluster_size_test/agpt/ncl_ntr/nfitstruc_32/nfits_10/n2body_128/grow_1.0/run_10_15/'
file = 'prediction_errors.out'
j1file = open(dir+file,'r')
nHold=250
def nstrip(list):
    '''Strips off /n'''
    import string
    list2 = []
    for string1 in list:   
        string2 = string1.strip("\n")
        list2.append(string2)
    return list2
lines = nstrip(j1file.readlines())[-nHold:-1] #get last fit of the Nfit.  

for i,line in enumerate(lines):
    print i, line.split()
    print line.split()[-3], line.split()[-2], line.split()[-1]