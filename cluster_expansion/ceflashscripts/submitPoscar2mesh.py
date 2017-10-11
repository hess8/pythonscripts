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
import voidWeightingRelax

def getVCmesh(dir,targetNmesh,meshtype,params) :

    lastDir = os.getcwd()   
#     meshc = dynamicPacking7.dynamicPack() #instance
    meshc = voidWeightingRelax.voidWeight() #instance
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
# os.chdir('/fslhome/bch/cluster_expansion/vcmesh/semiconductors/10OctRnosearchPerp/r0/Si_/bcc_6')
os.chdir('/fslhome/bch/cluster_expansion/vcmesh/semiconductors/sc_SiLP/r0/Si_/bcc_5')
# os.chdir('/fslhome/bch/cluster_expansion/vcmesh/test/r0/Al_1/bcc_5')
# os.chdir('/fslhome/bch/cluster_expansion/vcmesh/semiconductors/28SepRelaxVoidsFine/r2/Si_/bcc_13')
 
ntarget = 5
type = 'bcc'

paramLabels = ['wallClose','useVoids',    'rcutoff',  'tooClose','tooPlanar' 'rvCutoff','vwPower','wallPower','relax','interPower','wallFactor','wallOffset']
params =         ['0.05',       '0',         '3.0',       '-1',      '-1',       '4.0'   ,  '3.0' ,     '2.0' ,     '1',       '5.0',    '1.3',      '0.0']    #0.50   3.00   1.00   0.500.   1.   3.  -1.  -1.   3.   3.   3.   1.   6.   0.5  0.

'''Best Si from master branch:  
params0 =     [ 5.0 ] 
params1 =     [ 2.0 ]#wallPower
params2 =     [ 1.3 ] #wallfactor   
params3 =     [ 0.05] #wallClose    
params4 =     [ 0.0 ] #wallOffset    
params5 =     [ 0.5 ] #dw'''

statusOK,nops = getVCmesh(os.getcwd(),ntarget,type,params) 
if statusOK:
    writefile([],'OK')  
elif os.path.exists('OK'):
    os.system('rm OK')
    
    
    