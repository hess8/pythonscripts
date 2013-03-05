import sys,os,subprocess,time,sys,shutil, numpy as np
#from ceScriptTools import runningJobs
dir = '/fslhome/bch/cluster_expansion/cluster_size_test/agpt/ncl_ntr/nfitstruc_32/nfits_10/n2body_128/grow_1.0/run_10_15/'
file = 'J.1.out'
#file = 'clusters.out'
nclust = 200
clOrder = np.zeros((nclust),dtype=int)
clVert = np.zeros((nclust,6,3), dtype=float) #vertices 

j1file = open(dir+file,'r')
lines = j1file.readlines()
j1file.close()
clIndex = 0
for i,line in enumerate(lines):
#    print i, line
    if 'Number of vertices' in line:       
        print i, lines[i+1]
        order = int(lines[i+1])
        clOrder[clIndex] = order
        print "order", order
        for j in range(order):
            [x,y,z] = lines[i+7+j].split()[:3] #x, y, z coordinates
            clVert[clIndex,j,0]=x
            clVert[clIndex,j,1]=y
            clVert[clIndex,j,2]=z
        clIndex += 1           
print clIndex
print clOrder
    

