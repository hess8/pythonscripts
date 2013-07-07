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

def writejobname(jobfile,jobname):
    job = open(jobfile,'r')
    lines1 = job.readlines()
    job.close()
    for i,line in enumerate(lines1):
        if 'job-name' in line:
            lines1[i]='#SBATCH --job-name=%s\n' % jobname
    job = open(jobfile,'w')
    job.writelines(lines1)
    job.close()

maxDays = 4 #days to run this script
waitMin = 0.2 #minutes between checking 
nQueueWanted = 0 #keep this many in queue
#jobname = '3Gb4p' #best to update this
#jobname = '1Gb10proc'
user = 'bch'
mainDir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA4Gb4proc2/'
os.chdir(mainDir)
#mainDir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA1GB10proc/'
jobfile = mainDir + 'aflowjob'
#writejobname(jobfile,jobname)
toRunFile = mainDir + 'jobs2run'
starttime = time.time()
 
#os.system('rm slurm-*.out')     
while time.time()-starttime < maxDays*3600*24: 
    pending = subprocess.check_output(['squeue','-u',user,'--state=PENDING'])
    running = subprocess.check_output(['squeue','-u',user,'--state=RUNNING'])    
#    pending = subprocess.check_output(['squeue','-u',user,'-n',jobname,'--state=PENDING'])
#    running = subprocess.check_output(['squeue','-u',user,'-n',jobname,'--state=RUNNING'])
#    locked = subprocess.check_output(['find','-name','LOCK','|','wc','-l']) #gives error
    print pending
    pending = pending.splitlines()
    running = running.splitlines() 
    if jobsleft(toRunFile) < 1:
        break #no files to run 
    nsubmit = min([nQueueWanted, jobsleft(toRunFile)-len(running)-len(pending)+2])
    for ijob in range(nsubmit):
        subprocess.call(['sbatch', jobfile])
    print
    print time.asctime(time.localtime())
    print 'Jobs that were pending: ',len(pending)-1 #one line is header     
    print 'Jobs that were running: ',len(running)-1 #one line is header 
    print 'Just submitted %s jobs via %s' % (nsubmit, jobfile)
    print 'Days remaining for this script to run:', round(maxDays-(time.time()-starttime)/3600/24,4)
    print 'Files remaining in job list:' , jobsleft(toRunFile)
    print 'Will check again in %s min' % waitMin
#    print 'Locked files:' , locked   
    print
    time.sleep(waitMin*60) #seconds between checks
print 'done'