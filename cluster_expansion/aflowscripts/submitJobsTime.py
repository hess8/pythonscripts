#!/usr/bin/python
import time, os, subprocess

####### This "if" block is for python versions lower than 2.7. Needed for subprocess.check_output. 
if "check_output" not in dir( subprocess ): # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f
####### above needed for python versions lower than 2.7.
    
def jobsleft(toRunFile):
    '''Returns number of lines in file'''
    jobsfile = open(toRunFile,'r')
    lines1 = jobsfile.readlines()
    jobsfile.close()
    return len(lines1)

maxDays = 2 #days to run this script
waitMin = 1 #minutes between checking 
nQueueWanted = 1 #keep this many in queue
user = 'bch'
mainDir = '/fslhome/bch/cluster_expansion/alir/'
jobfile = mainDir + 'aflowjob'
toRunFile = mainDir + 'jobs2run'
starttime = time.time()
os.chdir(mainDir)    
while time.time()-starttime < maxDays*3600*24: 
    pending = subprocess.check_output(['squeue','-u',user, '--state=PENDING'])
    running  = subprocess.check_output(['squeue','-u',user, '--state=RUNNING'])
    print pending
    pending = pending.splitlines()
    running = running.splitlines() 
    if jobsleft(toRunFile) < 1:
        break #no files to run 
    for ijob in range(min([nQueueWanted - len(pending)+1,jobsleft(toRunFile)])):
        subprocess.call(['sbatch', jobfile])
    print 'Number of jobs that were pending: ',len(pending)-1 #one line is header     
    print 'Just submitted %s jobs via %s' % (str(nQueueWanted - len(pending)+1), jobfile)
    print 'Number of jobs running: ',len(running)-1 #one line is header 
    print 'Days remaining for this script to run:', round(maxDays-(time.time()-starttime)/3600/24,4)
    print 'Number of files that have not finished processing:' , jobsleft(toRunFile)
    print
    time.sleep(waitMin*60) #seconds between checks
print 'done'