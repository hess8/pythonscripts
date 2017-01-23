
import sys,os,subprocess,time
import numpy as np
from numpy import pi
#sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/')

#from kmeshroutines import checkq


def checkq(user):
    maxDays = 30 #days to run this script
    waitMin = 0.2 #minutes between checking         
    starttime = time.time()    
    while time.time()-starttime < maxDays*3600*24:
#        alias sq='squeue -o "%.7i %.9P %.18j %.8u %.2t %.10M %.6D %R" -u bch' 
        
        pending = subprocess.check_output(['squeue','-o','%.7i %.9P %.18j %.8u %.2t %.10M %.6D %R','-u',user,'--state=PENDING'])
        running = subprocess.check_output(['squeue','-o','%.7i %.7P %.18j %.8u %.2t %.10M %.6D %R','-u',user,'--state=RUNNING'])    
    #    locked = subprocess.check_output(['find','-name','LOCK','|','wc','-l']) #gives error
        print running
        pending = pending.splitlines()
        running = running.splitlines() 
        print 'Jobs that are pending: ',len(pending)-1 #one line is header     
        print 'Jobs that are running: ',len(running)-1 #one line is header 
        print 'Will check again in %s min' % waitMin
    #    print 'Locked files:' , locked   
        print
        time.sleep(waitMin*60) #seconds between checks

checkq('bch') #loops for days
