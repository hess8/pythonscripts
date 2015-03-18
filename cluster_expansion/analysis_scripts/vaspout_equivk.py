#!/usr/bin/python
'''
n is a multiple (m) of the volume factor  
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

#title_detail =  'Si:Si, -Accurate- precision'
#title_detail =  'Si:Si,no symmetry,cubic mesh,f1-50,ediff 1e-7 '
title_detail =  'Al:Al, cubic mesh '
#title_detail =  'Cu:Cu, cubic mesh,Accurate prec.  '
testfile = 'POSCAR'

def getibest(dirs):
#    mrange = []
    mmax = 0
    for i,dir in enumerate(dirs):
        if dir[1]=='1' and dir[2]=='_':
            m = int(dir.split('_')[-1])
            if m > mmax:
                ibest = i
                mmax = m
#        print 'CHOOSING F1_42 for test'
#        if dir == 'f1_42':
#            ibest = i
#            mmax = int(dir.split('_')[-1])
            
    return mmax, ibest

################# script #######################
#path = '/fslhome/bch/cluster_expansion/alal/cubic_al/equivk_c1-6_accurate/'
#path = '/fslhome/bch/cluster_expansion/alal/equivk/'
#path = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6_encut500/'
#path = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6.prec.accurate/'
#path = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6.tetra/'
#path =  '/fslhome/bch/cluster_expansion/alal/equivk.fixNGs.1_24.tetra/'
path =  '/fslhome/bch/cluster_expansion/alal/equivk_f-16.tetra.noBlochl/'
#path = '/fslhome/bch/cluster_expansion/cucu/equivk_accurate/'
#path ='/fslhome/bch/cluster_expansion/sisi/equivk_accurate/'
#path = '/fslhome/bch/cluster_expansion/alal/equivk_f1-6.acc.fixEF/'

#cubdir = path + 'structs.cubmesh/' #for plotting comparison
#for maindir in [path + 'structs.cubmesh/',path + 'structs.fccmesh/']:
for maindir in [path + 'structs.cubmesh/']:
#for maindir in [path + 'structs.fccmesh/']:
#path + 'structs.cubmesh/',
#path + 'structs.fccmesh/'
#]:

    #reallatt = zeros((3,3))
    os.chdir(maindir)
    dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
    if 'structselect' in dirs: dirs.remove('structselect')
    [mmax, ibest] = getibest(dirs)
    print 'energy of structure 1, multiplier %i, index %i used as ebest' % (mmax, ibest)

    #file1 = open('varypf.csv','a')
    #file1.write('Structure,Lattice,amax/amin,pfB,pf_orth,pf_orth2fcc,pf_maxpf, pf_pf2fcc, pfmax, meshtype' + ',' \
    #             + 'Improvement,fcc compatibility,Nmesh,TargetNmesh,Nmesh/Target,cbest' + '\n')
    #for i,directory in enumerate(dirs):    
    print dirs
    writeEnergiesOszicar(dirs) 
    writefermi(dirs) #skip so don't have to read it every time
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
    
    #try:
    #    file2 = open('lattype','r'); lattype = file2.read(); file2.close()
    #except:
    #    lattype = ''
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
#            print 'Run at location %s is not converged; removed' % i
    
    en_per_atom = array([float(energies[i])/dets[i]/2 for i in range(len(dets))]) #note: for silicon must divide by 2!!!!!!!!!!!!!!!!!!!!11
    ebest = en_per_atom[ibest]
    efermis = array([float(efermis[i]) for i in range(len(efermis))])
    efbest = float(efermis[ibest])
    print 'ebest', ebest
    print 'e-fermi best', efbest
    #for spreadsheet:
    err = abs(en_per_atom-ebest) 
    ef_err = abs(efermis - efbest) 
    ns = [ms[j]*dets[j] for j in range(len(dets))]  
    for i in range(len(names)):
        linei = names[i]+','+str(ns[i])+','+energies[i]+','+str(efermis[i])+','+str(dets[i])+','+str(en_per_atom[i])+','+str(err[i])+','+str(ef_err[i])+','+elconverge[i]+','+elsteps[i]+','+NkIBZ[i]+','+ str(round(float(cputime[i])/60,2))+'\n'        
        outfile.write(linei)    
    outfile.close() 
    
    struct = maindir.split('/')[-2]
    #plots
    [[ns,ms,en_per_atom,efermis,NkIBZ,NkfullBZ,dets],zerolist] = removezeros([ns,ms,en_per_atom,efermis,NkIBZ,NkfullBZ,dets])#for structures that were not finished and have zero entries



#    ebest = en_per_atom[argmax(ns)]
#    efbest = efermis[argmax(ns)]
    names2 = deepcopy(names); names = delete(names2,zerolist,axis=0)
    parts2 = deepcopy(parts); parts2 = delete(parts,zerolist,axis=0)
    parts3 = array([parts2[i]/dets[i] for i in range(len(dets))])
    
    ibestnew = nonzero(en_per_atom == ebest)[0][0]
    print 'New index for ebest', ibestnew
#    for plots
    err = abs(en_per_atom-ebest)+1e-9 #do these again with only the finished runs
    ef_err = abs(efermis - efbest)+1e-9
    #Write to file to simplify plotting for comparisons with other runs
    errlist = [str(ns[i])+'  '+str(err[i])+'\n' for i in range(len(ns))]
    print errlist
    writefile(errlist,'errlist_vs_n')

    #en_per_atom vs ns  
    titleadd = ''+ title_detail  
#    plotxy(ns,en_per_atom,'en_per_atom', titleadd + 'Vasp energy vs n (defines grid)','n','eV')
    
#    print ns
#    print err

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
    xlim((8.06, 8.08))  
    fig.savefig('ef_log_err_zoomed')  
    
    #log(ef_err) vs NkfullBZ 
    fig = figure()
    loglog(NkfullBZ,ef_err,'ro')
    title(titleadd + ' Error e-fermi')
    xlabel('Nk in unreduced BZ')
    ylabel('error in Ef') 
#    ylim((8.07, 8.08))  
    fig.savefig('err_ef_loglog_NkfullBZ')
    
        
    #log(ef_err) vs log(err)
    fig = figure()
    loglog(err,ef_err,'ro')
    title(titleadd + ' Error e-fermi')
    xlabel('error in energy (eV)')
    ylabel('error in Ef (eV)')   
    fig.savefig('err_ef_loglog_err')

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
    
        
    #log(err_EF) vs log(NkIBZ)
    fig = figure()
    loglog(NkIBZ,ef_err,'ro')
    title(titleadd + ' E-Fermi Error vs Nk in IBZKPT')
    xlabel('Nk')
    ylabel('error')   
    fig.savefig('nk_loglog_errEF') 
    
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
#    means = mean(parts3,axis=0)
#    deviations = abs(parts3 - means) + 1e-12 #last term to handle zero entries
    #print'deviations';print deviations
    
    partsbest = parts3[ibestnew,:]
    print 'partsbest',partsbest,'for the run', names[ibestnew]
#    sys.exit('stop')
#    partserr = zeros((len(dets),10))
    partsdev = zeros((len(dets),10))
    for j in range(9):
        for i in range(len(dets)):
            if abs(partsbest[j])>0:
#                partserr[i,j] = abs((parts3[i,j]-partsbest[j])/partsbest[j])
                partsdev[i,j] = abs(parts3[i,j]-partsbest[j])
            else:
#                partserr[i,j] = 1e-12
                partsdev[i,j] = 1e-12
    partsdev[:,9] = err

    fig = figure()
    ax1 = fig.add_subplot(111)
#   #Spreadsheet for error parts
#    outfile = open('parts_err_n.csv','w')
#    outfile.write('Struct,n_mesh,energy err,fermi err,alpha Z,Ewald energy,-1/2 Hartree,-exchange,-V(xc)+E(xc),PAW double counting,\
#                  entropy T*S,eigenvalues,atomic'+'\n')
#    for i in range(len(names)):
#        linei = names[i]+','+str(ns[i])+','+str(err[i])+','+str(ef_err[i])+','+str(partsdev[i,0])+','+str(partsdev[i,1])+','+str(partsdev[i,2])\
#        +','+ str(partsdev[i,3])+','+ str(partsdev[i,4])+','+ str(partsdev[i,5])+','+ str(partsdev[i,6])+','+ str(partsdev[i,7])+','+ str(partsdev[i,8])+'\n'        
#        outfile.write(linei)    
#    outfile.close() 
    
##        #Spreadsheet for error parts
#    outfile = open('parts_n.csv','w')
#    outfile.write('Struct,n_mesh,energy err,fermi err,alpha Z,Ewald energy,-1/2 Hartree,-exchange,-V(xc)+E(xc),PAW double counting,\
#                  entropy T*S,eigenvalues,atomic'+'\n')
#    for i in range(len(names)):
#        linei = names[i]+','+str(ns[i])+','+str(err[i])+','+str(ef_err[i])+','+str(parts3[i,0])+','+str(parts3[i,1])+','+str(parts3[i,2])\
#        +','+ str(parts3[i,3])+','+ str(parts3[i,4])+','+ str(parts3[i,5])+','+ str(parts3[i,6])+','+ str(parts3[i,7])+','+ str(parts3[i,8])+'\n'        
#        outfile.write(linei)    
#    outfile.close() 
#            
#    ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])

    enerlabels = ['alpha Z ','Ewald energy','-1/2 Hartree','-exchange','-V(xc)+E(xc)','PAW double counting',\
                  'entropy T*S','eigenvalues','atomic','total']
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
    #  atomic energy  EATOM  =       107.07526000 #not used!
    N = size(partsdev,axis = 1)    
    for i in range(N):
        ax1.semilogy(ns, partsdev[:,i],label=enerlabels[i],linestyle='None',color=cm.jet(1.*(i+1)/N), marker = 'o') # marker = 'o',
    plt.legend(loc='lower right',prop={'size':8});
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
        ax1.loglog(NkfullBZ, partsdev[:,i],label=enerlabels[i],linestyle='None',color=cm.jet(1.*(i+1)/N), marker = 'o') # marker = 'o',
    plt.legend(loc='upper right',prop={'size':6});
    show()
    fig.savefig('enerparts_err_NKfullBZ') 
    
#        plotxy(ns,parts3[:,7],'eigen_dev_vs_n', titleadd + 'Eigenvalue deviation vs n','n','eV')
  
#    eigenbest = parts3[-1,7] #take the last eigenvalue contribution to total energy as best
    
#    plotxy(ns,parts3[:,7],'eigen_dev_vs_n', titleadd + 'Eigenvalue dev vs n','n','eV')
#    fig = figure()
#    semilogy(ns,eigenerr,'ro')
#    title(titleadd + 'Eigenvalue err vs n')
#    xlabel('n')
#    ylabel('error')
#    fig.savefig('eigerr_log-n')  
    
    
    #plotxy(Nk,NkIBZ,'NkIBZ', titleadd+' Kpoint numbers','Nk from det M','Nk in IBZKPT')
#    plotxy(xaxis,cputime,'cpu', titleadd + ' CPU Time','Packing fraction','CPU time(sec)')
#    plotxy(NkIBZ,cputime,'cpu_vs_Nk', titleadd + ' CPU Time vs Nk in IBZKPT','Nk in IBZKPT','CPU time(sec)')

#*********** combine fcc with cubic mesh results *************

##title_detail = 'fcc mesh, Al:Al'
#testfile = 'POSCAR'
#os.chdir(cubdir)
#dirscub = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
#print dirscub
#writeEnergiesOszicar(dirscub)
#file = open('energies','r')
#energiescub = nstrip(file.readlines())
#file.close() 
#writedirnames(dirscub)
#detscub = getdata(dirscub,'detL')
#detscub = [float(detscub[i]) for i in range(len(detscub))]
#mscub = array(getms(dirscub),dtype = int)
#nscub = [str(float(mscub[i]*detscub[i])) for i in range(len(detscub))]
##print len(dirscub),len(energiescub), len(detscub)
#en_per_atomcub = array([float(energiescub[i])/detscub[i]/2 for i in range(len(detscub))])
#os.chdir(maindir)
#
##fcc and cubic comparison
#nscaled = [float(ms[i]*dets[i])*4**(1/3.0) for i in range(len(dets))] #scale fcc
#x1 = nscaled #scale fcc n to the cubic n 
#y1 = en_per_atom
#x2 = nscub
#y2 = en_per_atomcub
#from matplotlib.pyplot import *
#fig = figure()
#plot(x1, y1, 'ro',label = 'fcc mesh')
#plot(x2, y2, 'bo',label = 'cub mesh')
#title(titleadd + ' Vasp energy vs n (defines grid). n_fcc scaled to cubic')
#xlabel('n cubic, effective')
#ylabel('eV')
#legend(loc='lower right')
##ylim((-3.8,-3.7))
#show() 
#fig.savefig('vary_n_fcc_scaled_to_cub')  
#
##log plot of above
#errfcc = abs(en_per_atom-ebest) 
#errcub = abs(en_per_atomcub-ebest) 
##print errfcc
##print errcub
##print log10(errfcc)
##print log10(errcub)
#
#fig = figure()
#x1 = nscaled #scale fcc n to the cubic n 
#y1 = errfcc
#x2 = nscub
#y2 = errcub
##print 'nscaled', nscaled
##print 'errfcc', errfcc
#semilogy(x1, y1, 'ro',label = 'fcc mesh')
#semilogy(x2, y2, 'bo',label = 'cub mesh')
#title(titleadd + ' Vasp energy vs n (defines grid). n_fcc scaled to cubic')
#xlabel('n cubic, effective')
#ylabel('eV')
#legend(loc='lower right')
#ylim((1e-8,1e-2))
#xlim((0,40))
#show() 
#fig.savefig('vary_n_log_err_fcc_cub')  

#      
print 'Done'

