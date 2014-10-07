import os, subprocess, sys, time 


from matplotlib.pyplot import *
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/') 
from analysisToolsVasp import readfile 

file = 'structs.cubmesh/errlist_vs_n'
dirs = [
'/fslhome/bch/cluster_expansion/alal/equivk_f-16.tetra.noBlochl/', '/fslhome/bch/cluster_expansion/alal/equivk_f1-6.tetra/'
]
N = len(dirs)

#x = zeros[len(dirs)]
#y = zeros[len(dirs)]
labels = ['No correction','Blochl correction']
if len(labels) != len(dirs): sys.exit('Number of graphing labels is different from the number of datasets!\n  Stopping\n')

fig = figure()
ax1 = fig.add_subplot(111)
for idir,dir in enumerate(dirs):
    print idir
    lines = readfile(dir+file)
    x = [lines[i].strip().split()[0] for i in range(len(lines))]
    y = [lines[i].strip().split()[1] for i in range(len(lines))]
    ax1.semilogy(x, y,color=cm.jet(1.*(idir+1)/N),linestyle = 'None', marker = 'o', label=labels[idir])
#ylim((1e-12,1e0)) 
legend(loc='upper right',prop={'size':12});
title('Al:Al Tetrahedral methods\n Cubic mesh')
xlabel('n')
ylabel('error (eV)')
os.chdir(dirs[0]);os.chdir('../')
print dirs
fig.savefig('err_vs_n_multi')

print 'done'