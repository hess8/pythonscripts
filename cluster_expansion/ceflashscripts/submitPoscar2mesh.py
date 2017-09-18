#!/usr/bin/python

'''Define relaxed incommensurate kpoint mesh for this directory, and write INCAR'''
    
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,pi
from copy import copy, deepcopy
from numpy.linalg import norm
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/symmetry_k_mesh_search')
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar,create_poscar,readfile,writefile
# import dynamicPacking7
import voidWeighting

def getVCmesh(dir,targetNmesh,meshtype,params) :

    lastDir = os.getcwd()   
#     meshc = dynamicPacking7.dynamicPack() #instance
    meshc = voidWeighting.voidWeight() #instance
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
    statusOK,nops = meshc.pack(latticevecs,reciplatt,totatoms,aTypes,postype,transpose(positions),targetNmesh,meshtype,dir,params)
    os.chdir(lastDir)
    return statusOK,nops
#script:
os.chdir('/fslhome/bch/cluster_expansion/vcmesh/semiconductors/13SepFullWeights/bestRun/Si_/bcc_6')
ntarget = 6
type = 'bcc'
paramLabels = ['wallClose','rcutoff','tooClose','tooPlanar']
params =  ['0.5','3','1.0','0.5'] #0.50   3.00   1.00   0.50
statusOK,nops = getVCmesh(os.getcwd(),ntarget,type,params) 
if statusOK:
    writefile([],'OK')  
elif os.path.exists('OK'):
    os.system('rm OK')           