#!/usr/bin/python
'''    
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean
from numpy.linalg import norm
from analysisToolsVasp import writeEnergiesOszicar, writedirnames, nstrip, writeNk, writeNkIBZ, \
  writeElConverge, writeElSteps, writeCPUtime, enerparts, getdata, readfile, writefile, \
  getms, writefermi, removezeros
sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/analysis_scripts/plotting/') 
from plotTools import plotxy
from pylab import *
from copy import deepcopy
fprec=float64

testfile = 'POSCAR'

def getibest(dirs):
#    mrange = []
    mmax = 0
    for i,dir in enumerate(dirs):
        if dir[1]=='1' and dir[2]=='_':
            m = int(dir.split('_')[-1])
#            print dir, m
#            if m > mmax:
#                ibest = i
#                mmax = m
            if m == 24:
                ibest = i
                mmax = m
#            mrange.append(m)
    return mmax, ibest

def ebest_avg(ns):


################# script #######################
################# script #######################v
################# script #######################

paths = ['/fslhome/bch/cluster_expansion/structure_noise/sisi/equivk_accurate/',
         '/fslhome/bch/cluster_expansion/structure_noise/alal/equivk_f1-6.prec.accurate/'\
              , '/fslhome/bch/cluster_expansion/structure_noise/cucu/equivk_f1-6.prec.accurate/',
              '/fslhome/bch/cluster_expansion/structure_noise/gaga/equivk_accurate/',
              '/fslhome/bch/cluster_expansion/structure_noise/lili/equivk_accurate/',
              '/fslhome/bch/cluster_expansion/structure_noise/pdpd/equivk_accurate/',
              '/fslhome/bch/cluster_expansion/structure_noise/ptpt/equivk_accurate/']
nplots = len(paths)
titles = []
errsList = [[]]*nplots
nsList = [[]]*nplots
            
for iplot,path in enumerate(paths):
    titles.append(path.split('/')[-3][0].upper()+path.split('/')[-3][1].lower())
    for maindir in [path + 'structs.cubmesh/']:
        os.chdir(maindir)
        dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
        print dirs
        if 'structselect' in dirs: dirs.remove('structselect')
        [mmax, ibest] = getibest(dirs)
        print 'energy of structure 1, multiplier %i, index %i used as ebest' % (mmax, ibest)
        writeEnergiesOszicar(dirs) 
#        writefermi(dirs) #skip if don't want to read it every time
        writedirnames(dirs)
        writeNkIBZ(dirs)
        #writeNk(dirs)
        writeElConverge(dirs)
        writeElSteps(dirs)
        writeCPUtime(dirs)
        parts = enerparts(dirs)
        dets = getdata(dirs,'detL') #supercell volume in terms of unit cell volume
        dets = [float(dets[i]) for i in range(len(dets))]
        NkfullBZ = getdata(dirs,'detM') #Nk without symmetry reduction
        lattypes = getdata(dirs,'lattype')
        ms = array(getms(dirs),dtype = int)
        ################# summary #################
        outfile = open('vary_n.csv','w')
        outfile.write('Struct_m,n_mesh,energy,efermi,Natoms,energy/atom,energy err,fermi err,el converged,el steps,NIBZ,cputime(min\n')
        #os.chdir(mainDir)
        file = open(maindir+'names','r')
        names = nstrip(file.readlines())
        file.close()
        
        file = open(maindir+'energies','r')
        energies = nstrip(file.readlines())
        file.close()
        
        file = open(maindir+'efermi','r')
        efermis = nstrip(file.readlines())
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
          
        for i,stri in enumerate(elconverge): #make sure nonconverged runs have 0 energy
            if stri == 'N':
                energies[i] = '0.0'       
        en_per_atom = array([float(energies[i])/dets[i] for i in range(len(dets))])
        ebest = en_per_atom[ibest]
        efermis = array([float(efermis[i]) for i in range(len(efermis))])
        efbest = float(efermis[ibest])
        print 'ebest', ebest
        print 'e-fermi best', efbest
        err = abs(en_per_atom-ebest)/abs(ebest)
        ef_err = abs(efermis - efbest)/abs(efbest)
        ns = [ms[j]*dets[j] for j in range(len(dets))]  
        for i in range(len(names)):
            linei = names[i]+','+str(ns[i])+','+energies[i]+','+str(efermis[i])+','+str(dets[i])+','+str(en_per_atom[i])+','+str(err[i])+','+str(ef_err[i])+','+elconverge[i]+','+elsteps[i]+','+NkIBZ[i]+','+ str(round(float(cputime[i])/60,2))+'\n'        
            outfile.write(linei)    
        outfile.close() 
        
        struct = maindir.split('/')[-2]
        #plots: must remove all that failed
        [[ns,ms,en_per_atom,efermis,NkIBZ,NkfullBZ,dets],zerolist] = removezeros([ns,ms,en_per_atom,efermis,NkIBZ,NkfullBZ,dets])#for structures that were not finished and have zero entries
        parts2 = deepcopy(parts); parts2 = delete(parts,zerolist,axis=0)
        parts3 = array([parts2[i]/dets[i] for i in range(len(dets))])
        errabs = abs(en_per_atom-ebest)
        err = abs(en_per_atom-ebest)/abs(ebest) #do these again with only the finished runs
        ef_err = abs(efermis - efbest)/abs(efbest) 
        errsList[iplot] = err
        nsList[iplot] = ns
        
        #en_per_atom vs ns  
        titleadd = ''+ titles[iplot]  
        plotxy(ns,en_per_atom,'en_per_atom', titleadd + 'Vasp energy vs n (defines grid)','n','eV')
        fig = figure()
        semilogy(ns,err,'ro')
        title(titleadd + ' Error vs n (defines grid)')
        xlabel('n')
        ylabel('error')
        fig.savefig('vary_n_log_err')  
        
        #log(err) vs efermi
        fig = figure()
        semilogy(efermis,err,'ro')
        title(titleadd + ' Error vs e-fermi')
        xlabel('e-fermi (ev)')
        ylabel('error')   
        fig.savefig('ef_log_err')  
    
        #log(err) vs efermi zoomed
        fig = figure()
        semilogy(efermis,err,'ro')
        title(titleadd + ' Error vs e-fermi')
        xlabel('e-fermi (ev)')
        ylabel('error') 
        xlim((8.07, 8.08))  
        fig.savefig('ef_log_err_zoomed')  
        
        #log(ef_err) vs NkfullBZ zoomed
        fig = figure()
        loglog(NkfullBZ,ef_err,'ro')
        title(titleadd + ' Error vs e-fermi')
        xlabel('Nk in unreduced BZ')
        ylabel('error in Ef') 
    #    ylim((8.07, 8.08))  
        fig.savefig('err_ef_loglog_NkfullBZ')
    
        #log(err) vs NkIBZ
        fig = figure()
        semilogy(NkIBZ,err,'ro')
        title(titleadd + ' Error vs Nk in IBZKPT')
        xlabel('Nk')
        ylabel('error')   
        fig.savefig('nk_log_err')  
        
        #log(err) vs log(NkIBZ)
        fig = figure()
        loglog(NkIBZ,err,'ro')
        title(titleadd + ' Error vs Nk in IBZKPT')
        xlabel('Nk')
        ylabel('error')   
        fig.savefig('nk_loglog_err') 
        
        #log(err) vs log(NkfullBZ)
        fig = figure()
        loglog(NkfullBZ,err,'ro')
        title(titleadd + ' Error vs Nk in unreduced BZ')
        xlabel('Nk in unreduced BZ')
        ylabel('error')   
        fig.savefig('NkfullBZ_loglog_err') 
        
        #eigenerr plot
        eigenbest = parts3[-1,7] #take the last eigenvalue contribution to total energy as best
        eigenerr = plotxy(ns,parts3[:,7],'eigen_dev_vs_n', titleadd + 'Eigenvalue err vs n','n','eV')
       
    ##  enerparts deviations (eV) plots
        #find deviation from average for each part
        means = mean(parts3,axis=0)
        deviations = abs(parts3 - means) + 1e-12 #last term to handle zero entries
        
        partsbest = parts3[-1,:]
        partserr = zeros((len(dets),8))
        partsdev = zeros((len(dets),8))
        for j in range(8):
            for i in range(len(dets)):
                if abs(partsbest[j])>0:
                    partserr[i,j] = abs((parts3[i,j]-partsbest[j])/partsbest[j])
                    partsdev[i,j] = abs(parts3[i,j]-partsbest[j])
                else:
                    partserr[i,j] = 1e-12
                    partsdev[i,j] = 1e-12
    
        fig = figure()
        ax1 = fig.add_subplot(111)
    #    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
    
        enerlabels = ['alpha Z ','Ewald energy','-1/2 Hartree','-exchange','-V(xc)+E(xc)','PAW double counting',\
                      'entropy T*S','eigenvalues','atomic energy']
        xlabel('grid n')
        ylabel('error (eV)') 
        ylim((1e-12,1e0))
        
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
        N = size(partserr,axis = 1)    
        for i in range(N):
            ax1.semilogy(ns, deviations[:,i],label=enerlabels[i],linestyle='None',color=cm.jet(1.*(i+1)/N), marker = 'o') # marker = 'o',
        plt.legend(loc='upper right',prop={'size':6});
        show()
        fig.savefig('enerparts_deV') 
        
        ##### same plot, but relative error and vs NkfullBZ and loglog   
        fig = figure()
        #rcParams['axes.color_cycle']=['r','g']
        ax1 = fig.add_subplot(111)
    #    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
        xlabel('Nk in unreduced BZ')
        ylabel('relative error') 
        ylim((1e-12,1e0))    
        for i in range(N):
            ax1.loglog(NkfullBZ, deviations[:,i],label=enerlabels[i],linestyle='None',color=cm.jet(1.*(i+1)/N), marker = 'o') # marker = 'o',
        plt.legend(loc='upper right',prop={'size':6});
        show()
        fig.savefig('enerparts_err_NKfullBZ') 

os.chdir('/fslhome/bch/cluster_expansion/structure_noise')
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('n in cubic grid')
ylabel('Error (eV)') 
title('Structure noise\nReference energies: avg 4 highest n\'s')
xlim((0,55))
#ylim((1e-12,1e0))
for i in range(nplots):    
    ax1.semilogy(nsList[i], errsList[i],label=titles[i],linestyle='None',color=cm.jet(1.*(i+1)/float(nplots)), marker = 'o') # marker = 'o',
plt.legend(loc='upper right',prop={'size':14});
show()
string = ''
for i in range(nplots):
    string += '{}_'.format(titles[i])
fig.savefig('log_err_vs_n_{}accuratef1-6'.format(string)) 

     
print 'Done'

