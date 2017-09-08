#!/usr/bin/python

'''Define relaxed incommensurate kpoint mesh for this directory, and write INCAR'''
    
import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,pi
from copy import copy, deepcopy
from numpy.linalg import norm
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts/symmetry_k_mesh_search')
#import kmeshroutines as km
from kmeshroutines import nstrip, readposcar,create_poscar,readfile,writefile
import dynamicPacking7

def getVCmesh(dir,targetNmesh,meshtype,params) :

    lastDir = os.getcwd()   
    meshc = dynamicPacking7.dynamicPack() #instance
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
os.chdir('/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_Sigrid2Sep17/r0/Si_/bcc_8')

ntarget = 2
type = 'bcc'
#        ['power','wallPower-force','wallfactor','wallClose','wallOffset','dw' ]
params =    ['4.00',   '2.0 ',        '0.1',      '0.05',      '0.00',   '1.0']
statusOK,nops = getVCmesh(os.getcwd(),ntarget,type,params) 
if statusOK:
    writefile([],'OK')  
elif os.path.exists('OK'):
    os.system('rm OK')           