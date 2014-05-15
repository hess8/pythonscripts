#!/usr/bin/python
''' Reads structure numbers from training_set_structures.dat, 
    and creates directories with Aflow.in files. '''
    
import sys,os
################# functions #######################
def readTrainStruc(dir):
    os.chdir(dir)
    file1 = open('training_set_structures50.dat','r')
    lines = file1.readlines()
    file1.close()
    structs = [line.split()[1] for line in lines] #2nd column
    return structs

def aflowCreateJobs(structs,atomic,finalDir):
    nsubmit = 10 #number of structures to submit in each execution of aflow
    for i in range(len(structs)/nsubmit+1):
        if (i+1)*nsubmit < len(structs)-1:
            sublist = structs[i*nsubmit:(i+1)*nsubmit]
        else:
            sublist = structs[i*nsubmit:]
        strucstr = ''
        for struc in sublist:
            strucstr += 'f%s,'% struc
        try:
            commstr= './aconvasp --aflow_proto ' + strucstr+ ':%s' % atomic
            os.system(commstr)
        except:
            print 'Error in executing %s' % commstr
#    os.system('cp -r AFLOWDATA/* %s' % finalDir)
#    os.system('rm -r AFLOWDATA/')

def otherPrep():
    os.system("perl -pi -e 's/KPPRA=6000/KPPRA=10000/' */*/*/aflow.in")
    os.system("perl -pi -e 's/SYM=ON/SYM=OFF/' */*/*/aflow.in")
    os.system("find `pwd` -name 'aflow.in' > jobs2run")
    
################# script #######################

mainDir = '/fslhome/bch/cluster_expansion/alir/'
finalDir = mainDir
atomic = 'Al:Ir'
structs = readTrainStruc(mainDir)
aflowCreateJobs(structs,atomic,finalDir)
otherPrep()
print 'Done'