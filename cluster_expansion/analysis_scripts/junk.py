import sys,os,subprocess,time,sys,shutil, numpy as np
sys.path.append("/fslhome/bch/pythonscripts/cluster_expansion/run_scripts/")
from ceScriptTools import runningJobs

list1 = runningJobs() 
print list1
#list1=que.split('\n')[1:]
#list2=[]
#for i,item in enumerate(list1):
#    try:
#        list2.append(item.split()[2]) #3rd column
#    except:
#        print item, "can't be split"
#print list2