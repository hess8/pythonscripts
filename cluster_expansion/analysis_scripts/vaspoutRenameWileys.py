#!/usr/bin/python
'''
Comparison plot for different k mesh methods
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean,sort,argsort
from numpy.linalg import norm
from analysisToolsVasp import getEnergy, getNkIBZ,getNatoms, readfile, writefile,electronicConvergeFinish
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/ceflashscripts')
from symmetry import get_lattice_pointGroup, get_spaceGroup #these have vectors as ROWS
from kmeshroutines import readposcar

# from plotTools import plotxy
from pylab import figure,cm,plot,xlabel,ylabel,title,loglog,semilogy,legend,rcParams,style,rc
from copy import deepcopy

def areEqual(x,y,eps):
    return abs(x-y)<eps

testfile = 'POSCAR'

# paths = ['/fslhome/bch/cluster_expansion/mpmesh/cu.pt.ntest/12fstrEK_fcc'
#         ,'/fslhome/bch/cluster_expansion/vcmesh/the99']

paths = ['/fslhome/bch/cluster_expansion/vcmesh/the99sym_newMethod','/fslhome/bch/cluster_expansion/mpmesh/mpPure_MPbch']
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/semicond','/fslhome/bch/cluster_expansion/mpmesh/semicond']
# paths = ['/fslhome/bch/cluster_expansion/vcmesh/test','/fslhome/bch/cluster_expansion/mpmesh/semicond']
extpath = '/fslhome/bch/cluster_expansion/vcmesh/mueller_mp_data'

filter = 'Al' #string must be in dir name to be included
filter2 = 'Al_1'
summaryPath = paths[0]
# summaryPath = '/fslhome/bch/cluster_expansion/vcmesh/cu17Jul17/'
# summaryPath = paths[1]
iplot = 0
maxCalcs = 0
#count the number of plots:

#external data is of the form extpath/atom_method/struct.csv.  The csv has energies vs nK
os.chdir(extpath)
atoms_methods = sorted([d for d in os.listdir(extpath) if os.path.isdir(d) and filter in d])# os.chdir(extpath)
for atom_method in atoms_methods:
    atom = atom_method.split('_')[0]
    os.chdir(atom_method)
    os.system('rm -r .*lock*')
    structs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and filter in d])
    for struct in os.listdir(os.getcwd()):
        os.system('mv {} {}_{}'.format(struct,atom,struct))   
    os.chdir(extpath)


print 'Done'

