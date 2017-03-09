#!/usr/bin/env python
''' For each structure dir in maindir, add POSCAR using aconvasp.  dir can be empty  '''
    
import sys,os
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/AlIr/'
maindir = '/fslhome/bch/cluster_expansion/vcmesh/cu.pt.ntest/AFLOWDATAn/Cu_pvPt/'
os.chdir(maindir)
dirs= sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
print dirs
atomic = 'Cu:Pt'
for dir in dirs:
    print dir
    os.chdir(dir)     
    commstr2= 'aconvasp --proto=%s:%s | aconvasp --poscar>POSCAR' % (dir,atomic)   #for POSCAR creation         
    os.system(commstr2)
    os.chdir(maindir)

################## functions #######################
#def readTrainStruc(dir,filename):
#    os.chdir(dir)
#    file1 = open(filename,'r')
#    lines = file1.readlines()
#    file1.close() 
#    structs = [line.split()[1] for line in lines] #2nd column
#    return structs
#
#def aflowCreateJobs(structs,atomic,finalDir):
#    nsubmit = 1 #number of structures to submit in each execution of aflow
#    for i in range(len(structs)/nsubmit+1):
#        if (i+1)*nsubmit < len(structs)-1:
#            sublist = structs[i*nsubmit:(i+1)*nsubmit]
#        else:
#            sublist = structs[i*nsubmit:]
#        strucstr = ''
#        for struc in sublist:
##            strucstr += 'f%s,'% struc 
#            print struc  
#            strucstr += 'f%s'% struc #note the comma, and long string of structures in 'f%s,' doesn't work for '--proto='    
##            strucstr += 'h%s,'% struc
#            
#        try:
##            'aconvasp --proto=f1:Ag:Zr | aconvasp --poscar>POSCAR'
##            commstr= 'aconvasp --aflow_proto ' + strucstr+ ':%s' % atomic 
#            commstr2= 'aconvasp --proto=%s:%s | aconvasp --poscar>POSCAR' % (strucstr,atomic)   #for POSCAR creation         
##            commstr= 'aconvasp --aflow_proto ' + strucstr+ ':%s' % atomic
##            os.system(commstr)
#            print commstr2
#
#            os.system(commstr2)
#        
#        except:
#            print 'Error in executing %s' % commstr
#    os.system('cp -r AFLOWDATA/* %s' % finalDir)
#    os.system('rm -r AFLOWDATA/')
#
#def otherPrep():
##    os.system("perl -pi -e 's/KPPRA=6000/KPPRA=10000/' */*/*/aflow.in")
##    os.system("perl -pi -e 's/SYM=ON/SYM=OFF/' */*/*/aflow.in")
#    os.system("find `pwd` -name 'aflow.in' > jobs2run")
#    
################## script #######################
##filename='training_set_structures.dat'
##filename='f11000.dat'
#filename='f1_50.dat'
#mainDir = '/fslhome/bch/cluster_expansion/alir/'
##mainDir = '/fslhome/bch/cluster_expansion/hexagonal/'
##finalDir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/'
##finalDir = '/fslhome/bch/cluster_expansion/alir/h1_50/'
#finalDir = '/fslhome/bch/cluster_expansion/alir/f1_50b/'
#if not os.path.isdir(finalDir):
#    os.system('mkdir %s' % finalDir)
#atomic = 'Al:Ir'
#structs = readTrainStruc(mainDir,filename)
#print structs
#aflowCreateJobs(structs,atomic,finalDir)
#os.chdir(finalDir)
#otherPrep()
print 'Done'