def vasp_finish_check_spreadsheet(maindir,run,enerIsolatedAd,nsites):

    import os,subprocess,math,time 
    from numpy import array, zeros, binary_repr,log2,round
    import numpy
    from analysisTools import *
    #maindir = '/fslhome/bch/cluster_expansion/hexagonal/2x2adatoms/'
    os.chdir(maindir)
    #run = 'relaxfinal'
    #run = 'relax'
    #enerIsolatedAd = -.45365181e+01 # for tungsten (W) only
    #nsites = 8 
    
    #standard energies:
    eIsolatedH = -1.115
    eIsolatedC = -1.3179
    energyGraphane = -25.63
    nPrimCells = nsites/2
    bindEnergyGraphane = nPrimCells * (energyGraphane - 2*eIsolatedH - 2* eIsolatedC)
    
    dirlist = allFoldersList(maindir,run)
    nstruct = len(dirlist)
    #initialize arrays
    ener = zeros(nstruct); distances = zeros(nstruct); binde = zeros(nstruct);binde_nAd = zeros(nstruct); avgdCC =zeros(nstruct); 
    avgzC = zeros(nstruct); avgdAdC = zeros(nstruct); avgzAd = zeros(nstruct); min_dxyAdC = zeros(nstruct); 
    binaryatoms = zeros(nstruct,dtype=numpy.int); nAd = zeros(nstruct,dtype=numpy.int);steps = zeros(nstruct,dtype=numpy.int)
    finish = zeros(nstruct, dtype=numpy.str)
    
    structs = []
    
    NSW = getNSW('../vaspinput/%s/' % run) 
    for i,dir in enumerate(dirlist):
        print dir
        os.chdir(dir)
        struct = dir.split('/')[-3] ; structs.append(struct)
        ener[i] = energyOszicar()
        binaryatoms[i] = bin(int(struct.replace('struct','')))[2:]  #number only
        print binaryatoms[i]
        nAd[i] = nadatoms(str(binaryatoms[i]))
        nH = nsites-nAd[i]
        benergy_system = ener[i] -  nAd[i] * enerIsolatedAd - nH*eIsolatedH - nsites*eIsolatedC 
        binde[i] = benergy_system - bindEnergyGraphane #relative to graphane 
        if nAd[i] > 0:
            binde_nAd[i] = round(binde[i]/nAd[i],2)
        BEString = 'BE.vs.graphane'
        if FinishCheck():
            finish[i] = 'Y'
        else:
            finish[i] ='N'
        steps[i] = getSteps()
        if FinishCheck() and steps[i] < NSW:
            os.system('date > converged.dat')
        if os.path.exists('CONTCAR') and os.path.getsize('CONTCAR') > 0:       
            [avgdCC[i], avgzC[i], avgdAdC[i], avgzAd[i], min_dxyAdC[i]] = getPositionInfo()
        else:
            print 'CONTCAR not found, or empty'
    
    os.chdir(maindir)
    outfile = open('2x2all_%s.csv' %run,'w')
    outfile.write('Structure,Binary,%s,BE/nAdatoms,Calculated_Energy,N_Adatoms,AvgCCdist,AvgC_zcomp,AvgAdCdist,AvgAd_zcomp,OffTop_Dist,Finished,Steps\n' % BEString)
    
    #for i in range(len(elements)):
    
    #outfile.close()
    for i,dir in enumerate(dirlist):
        linei = structs[i]+','+ str(binaryatoms[i])+','+ str(binde[i]) +','+ str(binde_nAd[i]) +','+ str(ener[i])+','+ str(nAd[i])+','  \
            +str(avgdCC[i])+','+str(avgzC[i])+','+str(avgdAdC[i])+','+str(avgzAd[i])+','+str(min_dxyAdC[i])+ ','+str(finish[i])+','+ str(steps[i])+'\n'
        outfile.write(linei)
    
    outfile.close()