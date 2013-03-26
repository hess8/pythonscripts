''' Reads structure numbers from training_set_structures.dat, 
    and creates directories with Aflow.in files'''
    
import sys,os,subprocess
################# functions #######################
def readTrainStruc(dir):
    os.chdir(dir)
    file1 = open('training_set_structures.dat','r')
    lines = file1.readlines()
    file1.close()
#    print lines
    structs = [line.split()[1] for line in lines] #2nd column
    return structs

def aflowCreateJobs(structs,atomic,finalDir):
    for struct in structs:
        try:
            commstr= './aconvasp --aflow_proto f%s %s' % (struct, atomic)
            os.system(commstr)
        except:
            print 'Error in executing %s' % commstr
#    os.system('mv AFLOWDATA/* %s' % finalDir)
    
    
################# script #######################

mainDir = '/fslhome/bch/cluster_expansion/aflow1/test/'
finalDir = '/fslhome/bch/cluster_expansion/aflow1/stefano/LIBS/LIBtestbch/'
atomic = 'Al Pd,Pt'
structs = readTrainStruc(mainDir)
aflowCreateJobs(structs,atomic,finalDir)
print 'Done'


