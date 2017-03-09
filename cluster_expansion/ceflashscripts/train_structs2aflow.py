#!/usr/bin/env python
''' Reads structure numbers from training_set_structures.dat,
and creates directories with Aflow.in files. '''
    
import sys,os
################# functions #######################
def readTrainStruc(dir,filename):
    os.chdir(dir)
    file1 = open(filename,'r')
    lines = file1.readlines()
    file1.close()
    structs = [line.split()[1] for line in lines] #2nd column
    return structs

def aflowCreateJobs(structs,atomic,finalDir):
    nsubmit = 50 #number of structures to submit in each execution of aflow
    for i in range(len(structs)/nsubmit+1):
        if (i+1)*nsubmit < len(structs)-1:
            sublist = structs[i*nsubmit:(i+1)*nsubmit]
        else:
            sublist = structs[i*nsubmit:]
        strucstr = ''
        for struc in sublist:
# strucstr += 'f%s,'% struc
            strucstr += 'f%s,'% struc
        try:
            commstr= 'aconvasp --aflow_proto ' + strucstr+ ':%s' % atomic
            os.system(commstr)
        except:
            print 'Error in executing %s' % commstr
    os.system('cp -r AFLOWDATA/* %s' % finalDir)
    os.system('rm -r AFLOWDATA/')

def otherPrep():
# os.system("perl -pi -e 's/KPPRA=6000/KPPRA=10000/' */*/*/aflow.in")
# os.system("perl -pi -e 's/SYM=ON/SYM=OFF/' */*/*/aflow.in")
    os.system("find `pwd` -name 'aflow.in' > jobs2run")
    
################# script #######################
#filename='training_set_structures.dat'
#filename='f11000.dat'
filename='training_set_structures2.dat'
mainDir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/'
#finalDir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/'
finalDir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/AFLOWDATAn/'
if not os.path.isdir(finalDir):
    os.system('mkdir %s' % finalDir)
atomic = 'Cu:Pt'
#atomic = 'Al:Al'
structs = readTrainStruc(mainDir,filename)
print structs
aflowCreateJobs(structs,atomic,finalDir)
otherPrep()
#os.chdir(mainDir)
#print 'Done'