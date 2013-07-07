'''Parameters: processors, memory/proc
Script creates folder,jobname, edits jobfile, creates jobs2run, copies directory of structures.
Checks only type o '''

mkdir AFLOWDATA1Gb16proc/
cd AFLOWDATA1Gb16proc/

find `pwd` -name aflow.in>jobs2run
    pending = subprocess.check_output(['squeue','-u',user,'-n',jobname,'--state=PENDING'])
    running = subprocess.check_output(['squeue','-u',user,'-n',jobname,'--state=RUNNING'])