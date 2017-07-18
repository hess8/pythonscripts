import os,subprocess,time
'''Kills jobs in the 'jobrange' '''
user = 'bch'
# 
# jobrange = [   0, 1e10]  #for use with jobstring criterium
# jobstring = 'iso'  
# state = 'running'
# state = 'pending'
state = 'all'
jobrange = [ 16684877 , 16684900]  #must adjust this range of job numbers to kill
jobstring = ':'
joblist = subprocess.check_output(['squeue','-o','%.25i %.7P %.35j %.5u %.2t %.10M %.6D %R','-u',user,'--state={}'.format(state)])#'--state=pending'])
# joblist = subprocess.check_output(['squeue','-o','%.10i %.7P %.15j %.5u %.2t %.10M %.6D %R','-u',user])#'--state=pending']) 
#joblist = subprocess.check_output(['squeue','-o','%.7i %.7P %.15j %.5u %.2t %.10M %.6D %R','-u',user])#'--state=RUNNING'])    
##print joblist.split('\n')[1:][0].split()[0]
deleteIDs = []
for line in joblist.split('\n')[1:-1]: #the last and first lines should be ignored
#    
#    print int(line.split()[0])
    job = int(line.split()[0])
    name = line.split()[2]
#    if jobrange[0] <= job <= jobrange[1]:
    if jobrange[0] <= job <= jobrange[1] and jobstring in name:
        print line
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