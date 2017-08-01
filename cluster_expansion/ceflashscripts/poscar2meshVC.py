#!/usr/bin/python

'''Define relaxed incommensurate kpoint mesh for this directory, and write INCAR'''
    
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint
from copy import copy, deepcopy
from numpy.linalg import norm
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/symmetry_k_mesh_search')
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar,create_poscar,readfile,writefile
import dynamicPacking6



# 

def getVCmesh(dir,method,targetNmesh,meshtype) :

    lastDir = os.getcwd()   
    meshc = dynamicPacking6.dynamicPack() #instance
    [descriptor, scale, latticevecs, reciplatt, natoms, postype, positions] = readposcar('POSCAR',dir)
#         create_poscar('POSCAR',descriptor, scale, latticevecs, natoms, postype, positions, path) #just to remove the scale problem
    os.chdir(dir)
    totatoms = sum(natoms)
    atype = 1
    aTypes = []
    for natom in natoms:
        for i in range(natom):
            aTypes.append(atype)
        atype += 1
    aTypes = array(aTypes)
    statusOK,nops = meshc.pack(latticevecs,reciplatt,totatoms,aTypes,postype,transpose(positions),targetNmesh,meshtype,dir,method)
    os.chdir(lastDir)
    return statusOK,nops

#script:
ntarget = int(sys.argv[1])
type = sys.argv[2]
method = 0 #placeholder
statusOK,nops = getVCmesh(os.getcwd(),method,ntarget,type) 
if statusOK:
    writefile([],'OK')  
elif os.path.exists('OK'):
    print'Status not OK'
    os.system('rm OK')           