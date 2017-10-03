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
# os.chdir('/fslhome/bch/cluster_expansion/vcmesh/semiconductors/19SepFullWeights/bestRun/Si_/bcc_6')
# os.chdir('/fslhome/bch/cluster_expansion/vcmesh/test/r0/Si_/bcc_5')
os.chdir('/fslhome/bch/cluster_expansion/vcmesh/test/r0/Al_1/bcc_5')
# os.chdir('/fslhome/bch/cluster_expansion/vcmesh/semiconductors/28SepRelaxVoidsFine/r2/Si_/bcc_13')

ntarget = 5
type = 'bcc'
<<<<<<< HEAD
paramLabels = ['wallClose','useVoids',    'rcutoff',  'tooClose','tooPlanar' 'rvCutoff','vwPower','wallPower','relax','interPower','wallFactor','wallOffset']
params =         ['0.5',       '1',         '3.0',       '-1',      '-1',       '4.0'   ,  '2.0' ,     '3.0' ,     '1',       '6.0',       '0.5',       '0.0']    #0.50   3.00   1.00   0.50

'''    params0 =     [ 0.5 ]   #wallClose
    params1 =     [ 1 ]   #useVoids
    params2 =     [ 3.0 ]   #rcutoff
    params3 =     [ -1 ]   #tooClose
    params4 =     [ -1 ]  #tooPlanar
    params5 =     [ 4.0  ]   #rvCutoff
    params6 =     [ 2.0]  #vwPower
    params7 =     [ 3.0 ]   #wallPower
    params8 =     [ 1 ]   #relax (boolean)
    params9 =     [ 6.0 ]  #interPower
    params10 =    [ 0.5 ]  #wallFactor
    params11 =    [ 0.0 ]  #wallOffset'''

#         1.00   1.00   3.00  -1.00  -1.00   3.00   3.00   3.00   1.00   6.00   0.50   0.00 > out
=======
paramLabels = ['wallClose','useVoids','rcutoff','tooClose','tooPlanar','NvoidClosePoints','vwPower','wallPower',  'relax' 'interPower', 'wallFactor' , 'wallOffset']
params =         ['0.5',       '1',   '3',       '1.0',      '0.5',       '8'      ,        '2.0' ,     '6.0' ,        '1',       '3.0',       '0.5',       '0.0']    #0.50   3.00   1.00   0.50
#  ['0.5',       '3',     '1.0',     '0.5',       '8'      ,    '1.5'] 
>>>>>>> 3f6d1561040b799f220851f582a60dc346d1c84c
statusOK,nops = getVCmesh(os.getcwd(),ntarget,type,params) 
if statusOK:
    writefile([],'OK')  
elif os.path.exists('OK'):
    os.system('rm OK')
    
    
    