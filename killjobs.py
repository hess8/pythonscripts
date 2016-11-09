import os,subprocess,time
'''Kills jobs in the 'jobrange' '''
user = 'bch'
jobrange = [ 8880451,8880483]  #must adjust this range of job numbers to kill
running = subprocess.check_output(['squeue','-o','%.7i %.7P %.15j %.5u %.2t %.10M %.6D %R','-u',user])#'--state=RUNNING']) 
#running = subprocess.check_output(['squeue','-o','%.7i %.7P %.15j %.5u %.2t %.10M %.6D %R','-u',user])#'--state=RUNNING'])    
##print running.split('\n')[1:][0].split()[0]
deleteIDs = []
for line in running.split('\n')[1:-1]: #the last and first lines should be ignored
    print line
    print int(line.split(' ')[0])
    job = int(line.split(' ')[0])
    if jobrange[0] <= job <= jobrange[1]:
        deleteIDs.append(job)

print deleteIDs
print 'KILL JOBS!'
print 'Are you OK with deleting the above list of {} jobs?'.format(len(deleteIDs))
ans = 'no'
while ans != 'delete':
    ans = raw_input("If so, type 'delete' ")

for job in deleteIDs:
#    subprocess.call(['scancel','-u',user,str(job)])
    subprocess.call(['scancel',str(job)])
print "done"   