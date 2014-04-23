#!/usr/bin/python
'''    
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy,vasputil_dosplot
from pylab import *
fprec=float64

################# script #######################
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test101x/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA500/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/test/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATA11000/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test2/f6720/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test10^2/f5172/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test.10xNk/f3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test.noshift/f3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test/f3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test10^3/f3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/test10^5/f3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/testSi/f3/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/testSi/silicon/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/testSi/CCf3/'
maindir = '/fslhome/bch/cluster_expansion/sisi/test10^3/sidet2/'
#maindir = '/fslhome/bch/cluster_expansion/sisi/test10^5/'

title_detail = ''
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr/'
#maindir = '/fslhome/bch/cluster_expansion/alir/AFLOWDATAf1_50e/AlIr34-50/'

testfile = 'POSCAR'

#reallatt = zeros((3,3))
os.chdir(maindir)
dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
#file1 = open('varypf.csv','a')
#file1.write('Structure,Lattice,amax/amin,pfB,pf_orth,pf_orth2fcc,pf_maxpf, pf_pf2fcc, pfmax, meshtype' + ',' \
#             + 'Improvement,fcc compatibility,Nmesh,TargetNmesh,Nmesh/Target,cbest' + '\n')
#for i,directory in enumerate(dirs):    
print dirs
writeEnergiesOszicar(dirs) 
writedirnames(dirs)
writeNkIBZ(dirs)
writeNk(dirs)
writeElConverge(dirs)
writeElSteps(dirs)
writeCPUtime(dirs)
parts = enerparts(dirs)


print 'Plotting DOS'
for d in dirs:
    path = maindir + d +'/'
    try:
        vasputil_dosplot([], ["DOSCAR"], path) #options, args, dir
        plt.close()
    except:
        print 'Fail:'+ path
    os.chdir(maindir)

#print'parts'; print parts
#find deviation from average for each part
means = mean(parts,axis=0)
#print'means';print means
deviations = abs(parts - means) + 1e-16 #last term to handle zero entries
#print'deviations';print deviations

try:
    file2 = open('lattype','r'); lattype = file2.read(); file2.close()
except:
    lattype = ''
################# summary #################
outfile = open('vary_pf.csv','w')
outfile.write('pf,energy,el converged, el steps, Nkpoints, NIBZ,cputime(min\n')
 

#os.chdir(mainDir)
file = open(maindir+'names','r')
names = nstrip(file.readlines())
file.close()

file = open(maindir+'energies','r')
energies = nstrip(file.readlines())
file.close()

file = open(maindir+'Nk','r')
Nk = nstrip(file.readlines())
file.close()

file = open(maindir+'NkIBZ','r')
NkIBZ = nstrip(file.readlines())
file.close()

file = open(maindir+'elsteps','r')
elsteps = nstrip(file.readlines())
file.close()
elsteps

file = open(maindir+'elconverge','r')
elconverge = nstrip(file.readlines())
file.close()

file = open(maindir+'cputime','r')
cputime = nstrip(file.readlines())
file.close()

for i in range(len(names)):
    linei = names[i]+','+energies[i] +','+ elconverge[i]+','+elsteps[i]+','+ Nk[i]+','+NkIBZ[i]+','+ str(round(float(cputime[i])/60,2))+'\n'        
    outfile.write(linei)    
outfile.close() 

struct = maindir.split('/')[-2]
#plots
titleadd = struct +','+ lattype +','+ title_detail 
plotxy(names,energies,'vary_pf', titleadd + 'Vasp energy vs packing fraction','Packing fraction','eV')
#plotxy(Nk,NkIBZ,'NkIBZ', titleadd+' Kpoint numbers','Nk from det M','Nk in IBZKPT')
plotxy(names,cputime,'cpu', titleadd + ' CPU Time','Packing fraction','CPU time(sec)')
plotxy(NkIBZ,cputime,'cpu_vs_Nk', titleadd + ' CPU Time vs Nk in IBZKPT','Nk in IBZKPT','CPU time(sec)')

#enerparts deviations plots
fig = figure()
#rcParams['axes.color_cycle']=['r','g']
ax1 = fig.add_subplot(111)
ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
N = size(deviations,axis = 1)
enerlabels = ['alpha Z ','Ewald energy','-1/2 Hartree','-exchange','-V(xc)+E(xc)','PAW double counting',\
              'entropy T*S','eigenvalues','atomic energy']
#
#  alpha Z        PSCENC =        -0.43400971
#  Ewald energy   TEWEN  =      -150.19868154
#  -1/2 Hartree   DENC   =        -0.23270923
#  -exchange  EXHF       =         0.00000000
#  -V(xc)+E(xc)   XCENC  =       -52.00285543
#  PAW double counting   =        71.87483513       -8.01524693
#  entropy T*S    EENTRO =        -0.00011188
#  eigenvalues    EBANDS =        24.48126738
#  atomic energy  EATOM  =       107.07526000

for i in range(N):
    ax1.semilogy(names, deviations[:,i],label=enerlabels[i],linestyle='None',color=cm.jet(1.*(i+1)/N), marker = 'o') # marker = 'o',
plt.legend(loc='lower right');
show()
fig.savefig('enerparts_dev') 
      
print 'Done'

