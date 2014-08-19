#!/usr/bin/python
'''  
Choose a few structures to plot convergence, rather than the entire folder  
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



def getibest(dirs):
#    mrange = []
    mmax = 0
    for i,dir in enumerate(dirs):
        if dir[1]=='1' and dir[2]=='_':
            m = int(dir.split('_')[-1])
            print dir, m
            if m > mmax:
                ibest = i
                mmax = m
#            mrange.append(m)
    return mmax, ibest



################# script #######################
################# script #######################
################# script #######################





#path = '/fslhome/bch/cluster_expansion/sisi/equivk/'
#path = '/fslhome/bch/cluster_expansion/sisi/nosymequivk/'
#path = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6_encut300/'
path = '/fslhome/bch/cluster_expansion/alal/cubic_al/equivk_c1-6_encut500/'
#path = '/fslhome/bch/cluster_expansion/cucu/equivk/'

cubdir = path + 'structs.cubmesh/' #for plotting comparison


#title_detail =  'Si:Si'
#title_detail =  'Si:Si,no symmetry,cubic mesh,f1-50,ediff 1e-7 '
title_detail =  'Cubic Al:Al (2.86 ang), cubic mesh, encut 500'
#title_detail =  'Cu:Cu, cubic mesh,f1-50,ediff 1e-7 '

structselect = ['c1','c3']#,'f7', 'f8','f9','f10']#need to have f1 in here


testfile = 'POSCAR'



for maindir in [path + 'structs.cubmesh/']:
#for maindir in [path + 'structs.cubtest/']:
    tempdir = maindir+'structselect/'


    #reallatt = zeros((
    os.chdir(maindir)
    dirs1 = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
    if not os.path.isdir(tempdir): os.makedirs(tempdir)
    os.chdir(tempdir) #put lists for each structure in here, writing over them as we go
    os.system('rm *')
    for structi in structselect: 
        dirs = []
        dirsfull = []
        for dir in dirs1:
            if structi+'_' in dir:
                dirs.append(dir)
                dirsfull.append(maindir+dir)
        print dirs
        if structi ==  structselect[0]:
            [mmax, ibest] = getibest(dirs)
            print 'energy of structure 1, multiplier %i, index %i used as ebest' % (mmax, ibest)  
        
        #ibest = 17 ###################TEMP  ONLY !!!!!!!!!!!!!!!!!!!!!
        writeEnergiesOszicar(dirsfull) 
        writefermi(dirsfull) #skip so don't have to read it every time
        writedirnames(dirsfull)
        writeNkIBZ(dirsfull)
        #writeNk(dirsfull)
        writeElConverge(dirsfull)
        writeElSteps(dirsfull)
        writeCPUtime(dirsfull)
        parts = enerparts(dirsfull)
        dets = getdata(dirsfull,'detL') #supercell volume in terms of unit cell volume
        dets = [float(dets[i]) for i in range(len(dets))]
        NkfullBZ = getdata(dirsfull,'detM') #Nk without symmetry reduction
        lattypes = getdata(dirsfull,'lattype')
        ms = array(getms(dirsfull),dtype = int)
        
        filesstructselect = sorted([f for f in os.listdir(os.getcwd())])

        #os.chdir(mainDir)
        
        file = open(tempdir+'energies','r')
        energies = nstrip(file.readlines())
        file.close()
        os.system('mv %s %s' % ('energies', 'energies'+structi))
        
        file = open(tempdir+'efermi','r')
        efermis = nstrip(file.readlines())
        file.close()
        os.system('mv %s %s' % ('efermi', 'efermi'+structi))
                   
        file = open(tempdir+'NkIBZ','r')
        NkIBZ = nstrip(file.readlines())
        file.close()
        os.system('mv %s %s' % ('NkIBZ', 'NkIBZ'+structi))        
        
        file = open(tempdir+'elsteps','r')
        elsteps = nstrip(file.readlines())
        file.close()
        elsteps
        os.system('mv %s %s' % ('elsteps', 'elsteps'+structi))
        
        file = open(tempdir+'elconverge','r')
        elconverge = nstrip(file.readlines())
        file.close()
        os.system('mv %s %s' % ('elconverge', 'elconverge'+structi))
        #keep from overwriting these files:

        for i,stri in enumerate(elconverge): #make sure nonconverged runs have 0 energy
            if stri == 'N':
                energies[i] = '0.0'
        #            print 'Run at location %s is not converged; removed' % i
        
        en_per_atom = array([float(energies[i])/dets[i] for i in range(len(dets))])
        efermis = array([float(efermis[i]) for i in range(len(efermis))])        
        if structi ==  structselect[0]: #comparing everything to it for now
            ebest = en_per_atom[ibest]
            efbest = float(efermis[ibest]) 
            print 'ebest', ebest
            print 'e-fermi best', efbest
        err = abs(en_per_atom-ebest)
        ef_err = abs(efermis - efbest)
        ns = [ms[j]*dets[j] for j in range(len(dets))]  
        
        #plots
        [[ns,ms,en_per_atom,efermis,NkIBZ,NkfullBZ,dets],zerolist] = removezeros([ns,ms,en_per_atom,efermis,NkIBZ,NkfullBZ,dets])#for structures that were not finished and have zero entries
        #    ebest = en_per_atom[argmax(ns)]
        #    efbest = efermis[argmax(ns)]
        err = abs(en_per_atom-ebest)+1e-8 #do these again with only the finished runs
        ef_err = abs(efermis - efbest)+1e-8

        #write important arrays to file
        writefile([str(j)+'\n' for j in ns],'ns%s' % structi)
        writefile([str(j)+'\n' for j in ms],'ms%s' % structi)
        writefile([str(j)+'\n' for j in NkIBZ],'NkIBZ%s' % structi)
        writefile([str(j)+'\n' for j in NkfullBZ],'NkfullBZ%s' % structi)
        writefile([str(j)+'\n' for j in dets],'dets%s' % structi)
        writefile([str(j)+'\n' for j in err],'err%s' % structi)
        writefile([str(j)+'\n' for j in ef_err],'ef_err%s' % structi)
        
# log err vs n for all structs
os.chdir(tempdir)
fig = figure()
ax1 = fig.add_subplot(111)
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
xlabel('n in cubic grid')
ylabel('Error (eV)') 
title('Structure noise '+title_detail+':\nTheoretical values: max n on struct 1; 1e-8 mark')
xlim((0,55))
#ylim((1e-12,1e0))
for i,structi in enumerate(structselect):  
    nlist = readfile('ns%s' % structi)
    errlist = readfile('err%s' % structi)
    ax1.semilogy(nlist, errlist,label=structi,linestyle='None',color=cm.jet(1.*(i+1)/len(structselect)), marker = 'o') # marker = 'o',
plt.legend(loc='upper right',prop={'size':14});
show()
fig.savefig('log_err_vs_n_structs') 


        
#    for structi in structselect:
#        for arr in [] 

        
#    
#    #en_per_atom vs ns  
#    titleadd = ''+ title_detail  
#    plotxy(ns,en_per_atom,'en_per_atom', titleadd + 'Vasp energy vs n (defines grid)','n','eV')        
#    
#    

#    
#    fig = figure()
#    semilogy(ns,err,'ro')
#    title(titleadd + ' Error vs n (defines grid)')
#    xlabel('n')
#    ylabel('error')
#    fig.savefig('vary_n_log_err')  
#    
#    #log(err) vs efermi
#    fig = figure()
#    semilogy(efermis,err,'ro')
#    title(titleadd + ' Error vs e-fermi')
#    xlabel('e-fermi (ev)')
#    ylabel('error')   
#    fig.savefig('ef_log_err')  
#    
#    #log(err) vs efermi zoomed
#    fig = figure()
#    semilogy(efermis,err,'ro')
#    title(titleadd + ' Error vs e-fermi')
#    xlabel('e-fermi (ev)')
#    ylabel('error') 
#    xlim((8.07, 8.08))  
#    fig.savefig('ef_log_err_zoomed')  
#    
#    #log(ef_err) vs NkfullBZ zoomed
#    fig = figure()
#    loglog(NkfullBZ,ef_err,'ro')
#    title(titleadd + ' Error vs e-fermi')
#    xlabel('Nk in unreduced BZ')
#    ylabel('error in Ef') 
#    #    ylim((8.07, 8.08))  
#    fig.savefig('err_ef_loglog_NkfullBZ')
#    
#    #log(err) vs NkIBZ
#    fig = figure()
#    semilogy(NkIBZ,err,'ro')
#    title(titleadd + ' Error vs Nk in IBZKPT')
#    xlabel('Nk')
#    ylabel('error')   
#    fig.savefig('nk_log_err')  
#    
#    #log(err) vs log(NkIBZ)
#    fig = figure()
#    loglog(NkIBZ,err,'ro')
#    title(titleadd + ' Error vs Nk in IBZKPT')
#    xlabel('Nk')
#    ylabel('error')   
#    fig.savefig('nk_loglog_err') 
#    
#    #log(err) vs log(NkfullBZ)
#    fig = figure()
#    loglog(NkfullBZ,err,'ro')
#    title(titleadd + ' Error vs Nk in unreduced BZ')
#    xlabel('Nk in unreduced BZ')
#    ylabel('error')   
#    fig.savefig('NkfullBZ_loglog_err') 
#    
#    #eigenerr plot
#    eigenbest = parts3[-1,7] #take the last eigenvalue contribution to total energy as best
#        eigenerr = plotxy(ns,parts3[:,7],'eigen_dev_vs_n', titleadd + 'Eigenvalue err vs n','n','eV')

print 'Done'

