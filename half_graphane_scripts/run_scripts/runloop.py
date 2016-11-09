import os,sys,time,subprocess
sys.path.append('/bluehome2/bch//pythonscripts/half_graphane_scripts/analysis_scripts/supercell_all_structures_analysis')
sys.path.append('/bluehome2/bch/pythonscripts/half_graphane_scripts/run_scripts/supercell_all_structures_run/')

import graphane_super_ad,vasp_output
from graphane_super_ad import submitsupercell
from vasp_output import vasp_finish_check_spreadsheet

maindir = '/fslhome/bch/cluster_expansion/hexagonal/2x2adatoms/'
run = 'relaxfinal'
user = 'bch'

Nsuper1 = 2 # N of supercells in two directions
Nsuper2 = 2
atomslist = ['C','H','W']
dAd = 2.2  # Adatom  distance from plane

maxDays = 4 #days to run this script
waitMin = 1 #minutes between checking 
loopHrs = 6 #hrs before restarting all jobs
scriptstart = time.time()
jobNs=[]    
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
while time.time()-scriptstart < maxDays*3600*24: #time to ru
    ####### Stop jobs in this dir #######
    # need to give dir name to job
    os.system('squeue -o "%.7i %.9P %.18j %.8u %.2t %.10M %.6D %R" -u bch')
    print 'Stopping jobs'
    os.system('scancel -u bch --state=RUNNING')
    
    ####### Check jobs #######
    nsites = 8 
    enerIsolatedAd = -.45365181e+01 # for tungsten (W) only
    print 'Analyzing runs'
    vasp_finish_check_spreadsheet(maindir,run,enerIsolatedAd,nsites)
    
    #######  submit jobs #######
    submitsupercell(maindir,run,atomslist,dAd,Nsuper1,Nsuper2)
        
    ####### show running jobs #######
    loopstart = time.time() 
    while time.time()-loopstart < loopHrs*3600: 
        pending = subprocess.check_output(['squeue','-u',user,'--state=PENDING'])
        running = subprocess.check_output(['squeue','-u',user,'--state=RUNNING'])   
        if len(pending)>85: #lenthg of line if empty pending queue
            print pending
        pending = pending.splitlines(); running = running.splitlines() 
        print
        print time.asctime(time.localtime())
        print 'Jobs pending: ',len(pending)-1 #one line is header     
        print 'Jobs running: ',len(running)-1 #one line is header 
        if len(jobNs) > 0: print 'Jobs vs time', jobNs      
        print 'Hrs before jobs are restarted', round(loopHrs-(time.time()-loopstart)/3600,2)
        print 'Days remaining for this script to run:', round(maxDays-(time.time()-scriptstart)/3600/24,4) 
        time.sleep(waitMin*60)
    jobNs.append(len(pending)+len(running)-2) #update the n