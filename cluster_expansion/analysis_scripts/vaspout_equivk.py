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
fprec=float64

#title_detail =  'Si:Si'
title_detail =  'Al:Al, cubic mesh,f3-50 '
testfile = 'POSCAR'

################# script #######################
#path = '/fslhome/bch/cluster_expansion/sisi/equivk/'
path = '/fslhome/bch/cluster_expansion/alir/enumtest/'
#path = '/fslhome/bch/cluster_expansion/alir/enumtest/test/'
cubdir = path + 'structs.cubmesh/' #for plotting comparison
for maindir in [path + 'structs.cubmesh/']:
#path + 'structs.cubmesh/',
#path + 'structs.fccmesh/'
#]:

    #reallatt = zeros((3,3))
    os.chdir(maindir)
    dirs = sorted([d for d in os.listdir(os.getcwd()) if os.path.isdir(d)])
    #file1 = open('varypf.csv','a')
    #file1.write('Structure,Lattice,amax/amin,pfB,pf_orth,pf_orth2fcc,pf_maxpf, pf_pf2fcc, pfmax, meshtype' + ',' \
    #             + 'Improvement,fcc compatibility,Nmesh,TargetNmesh,Nmesh/Target,cbest' + '\n')
    #for i,directory in enumerate(dirs):    
    print dirs
    writeEnergiesOszicar(dirs) 
#    writefermi(dirs) #skip so don't have to read it every time
    writedirnames(dirs)
    writeNkIBZ(dirs)
    #writeNk(dirs)
    writeElConverge(dirs)
    writeElSteps(dirs)
    writeCPUtime(dirs)
    parts = enerparts(dirs)
    dets = getdata(dirs,'detL')
    dets = [float(dets[i]) for i in range(len(dets))]
    lattypes = getdata(dirs,'lattype')
    ms = array(getms(dirs),dtype = int)
    
    #print'parts'; print parts
    #find deviation from average for each part
    means = mean(parts,axis=0)
    #print'means';print means
    deviations = abs(parts - means) + 1e-16 #last term to handle zero entries
    #print'deviations';print deviations
    
    #try:
    #    file2 = open('lattype','r'); lattype = file2.read(); file2.close()
    #except:
    #    lattype = ''
    ################# summary #################
    outfile = open('vary_n.csv','w')
    outfile.write('Struct,n_mesh,energy,efermi,Natoms, energy/atom,el converged, el steps, NIBZ,cputime(min\n')
    
    
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
    
    en_per_atom = array([float(energies[i])/dets[i] for i in range(len(dets))])
    ns = [str(float(ms[i]*dets[i])) for i in range(len(dets))]
    
    for i in range(len(names)):
        linei = names[i]+','+ns[i]+','+energies[i]+','+efermis[i]+','+str(dets[i])+','+str(en_per_atom[i])+','+elconverge[i]+','+elsteps[i]+','+NkIBZ[i]+','+ str(round(float(cputime[i])/60,2))+'\n'        
        outfile.write(linei)    
    outfile.close() 
    
    struct = maindir.split('/')[-2]

    #plots
    [ns,ms,en_per_atom,efermis,NkIBZ] = removezeros([ns,ms,en_per_atom,efermis,NkIBZ])#for structures that were not finished and have zero entries
    
    #en_per_atom vs ns  
    titleadd = ''+ title_detail  
    plotxy(ns,en_per_atom,'en_per_atom', titleadd + 'Vasp energy vs n (defines grid)','n','eV')
    
    #log(err) vs ns  
#    ebest = en_per_atom[argmax(ns)] ########## Need to find a better "best"
    N = len(en_per_atom); 
    ebest = average(en_per_atom)#[int(N/2):]) #take last half average as the best (largest n's)
    err = abs(en_per_atom-ebest)/abs(ebest)
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
    xlim((8.05, 8.1))  
    fig.savefig('ef_log_err_zoomed')  

    
    #log(err) vs NkIBZ
    fig = figure()
    semilogy(NkIBZ,err,'ro')
    title(titleadd + ' Error vs Nk in IBZKPT')
    xlabel('Nk')
    ylabel('error')   
    fig.savefig('nk_log_err')  
    
    #plotxy(Nk,NkIBZ,'NkIBZ', titleadd+' Kpoint numbers','Nk from det M','Nk in IBZKPT')
#    plotxy(xaxis,cputime,'cpu', titleadd + ' CPU Time','Packing fraction','CPU time(sec)')
#    plotxy(NkIBZ,cputime,'cpu_vs_Nk', titleadd + ' CPU Time vs Nk in IBZKPT','Nk in IBZKPT','CPU time(sec)')

##combine with cubic mesh results 
#
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
#en_per_atomcub = array([float(energiescub[i])/detscub[i] for i in range(len(detscub))])
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
#ebest = en_per_atomcub[argmax(nscub)]
##ebest = -10.8392895
#print 'ebest',ebest
#errfcc = abs(en_per_atom-ebest)/abs(ebest)
#errcub = abs(en_per_atomcub-ebest)/abs(ebest)
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
#semilogy(x1, y1, 'ro',label = 'fcc mesh')
#semilogy(x2, y2, 'bo',label = 'cub mesh')
#title(titleadd + ' Vasp energy vs n (defines grid). n_fcc scaled to cubic')
#xlabel('n cubic, effective')
#ylabel('eV')
#legend(loc='lower right')
#ylim((1e-8,1e-2))
##xlim((0,68))
#show() 
#fig.savefig('vary_n_log_err_fcc_cub')  
#
#
##enerparts deviations plots
#fig = figure()
##rcParams['axes.color_cycle']=['r','g']
#ax1 = fig.add_subplot(111)
#ax1.set_color_cycle(['r','b','g','c', 'm', 'y', 'k'])
#N = size(deviations,axis = 1)
#enerlabels = ['alpha Z ','Ewald energy','-1/2 Hartree','-exchange','-V(xc)+E(xc)','PAW double counting',\
#              'entropy T*S','eigenvalues','atomic energy']
##ylim((-3.7,-3.8))
#
##
##  alpha Z        PSCENC =        -0.43400971
##  Ewald energy   TEWEN  =      -150.19868154
##  -1/2 Hartree   DENC   =        -0.23270923
##  -exchange  EXHF       =         0.00000000
##  -V(xc)+E(xc)   XCENC  =       -52.00285543
##  PAW double counting   =        71.87483513       -8.01524693
##  entropy T*S    EENTRO =        -0.00011188
##  eigenvalues    EBANDS =        24.48126738
##  atomic energy  EATOM  =       107.07526000
#
#for i in range(N):
#    ax1.semilogy(ns, deviations[:,i],label=enerlabels[i],linestyle='None',color=cm.jet(1.*(i+1)/N), marker = 'o') # marker = 'o',
#plt.legend(loc='lower right');
#show()
#fig.savefig('enerparts_dev') 
#      
print 'Done'

