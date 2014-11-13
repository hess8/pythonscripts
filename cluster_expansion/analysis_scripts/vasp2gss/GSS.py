'''
Created on Aug 29, 2014

@author: eswens13
'''
import os, subprocess
from numpy import amin, amax, sort, zeros, array, floor, exp, ceil

class GSS:
    

    def __init__(self, atoms, volRange, plotTitle, xlabel, ylabel, uncleOutput):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.volRange = volRange
        self.plotTitle = plotTitle
        self.xlabel = xlabel
        self.ylabel = ylabel
        
        self.enumFolder = os.getcwd() + '/enum/'
        self.fitsDir = None
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
        self.Ncs = []
        self.priorities = []

    def readfile(self,filepath): #bch
        file1 = open(filepath,'r')
        lines = file1.readlines()
        file1.close()
        return lines
    
    def makeGSSDirectories(self):
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                self.fitsDir = atomDir + '/fits/'
                gssDir = atomDir + '/gss/'
                subprocess.call(['mkdir',gssDir])
                subprocess.call(['cp',self.enumFolder + 'struct_enum.out', gssDir])
                subprocess.call(['cp',self.enumFolder + 'lat.in', gssDir])
                subprocess.call(['cp',self.fitsDir + 'J.1.out', gssDir])
                
                infile = open(self.neededFilesDir + 'groundstatesearch.in','r')
                inlines = [line for line in infile]
                infile.close()
                
                outfile = open(gssDir + 'groundstatesearch.in','w')
                for i in xrange(len(inlines)):
                    if i == 4:
                        outfile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + "\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
                
                infile = open(self.neededFilesDir + 'gss_plot.gp','r')
                inlines = [line for line in infile]
                infile.close()
                
                outfile = open(gssDir + 'gss_plot.gp','w')
                for i in xrange(len(inlines)):
                    if i == 3:
                        outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atom) + "\"\n")#bch
                    elif i == 4:
                        outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                    elif i == 5:
                        outfile.write("set title \"" + self.plotTitle + " (" + atom + ")\"\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
                
                #Binding energy bch
                infile = open(self.neededFilesDir + 'BE_plot.gp','r')
                inlines = [line for line in infile]
                infile.close()
                outfile = open(gssDir + 'BE_plot.gp','w')
                for i in xrange(len(inlines)):
                    if i == 3:
                        outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atom) + "\"\n")#bch
                    elif i == 4:
                        outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                    elif i == 5:
                        outfile.write("set title \"" + 'Binding energy vs graphene'+ " (" + atom + ")\"\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()

                #Hexagonal formation energy energy bch        
                data = self.readfile(atomDir+'/vaspHFE.out')
                ydata = [float(line.strip().split()[1]) for line in data]
                ymax = min(amax(ydata),1.0) #don't let yrange get over 1 eV
                infile = open(self.neededFilesDir + 'HFE_plot.gp','r')
                inlines = [line for line in infile]
                infile.close()
                outfile = open(gssDir + 'HFE_plot.gp','w')
                for i in xrange(len(inlines)):
                    if i == 3:
                        outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atom) + "\"\n")#bch
                    elif i == 4:
                        outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                    elif i == 5:
                        outfile.write("set title \"" + 'Formation energy vs H2, metal hex monolayer'+ " (" + atom + ")\"\n")
                    elif 'plot "' in inlines[i]:         
                        outfile.write('set yrange [:{}]\n'.format(ymax))
                        outfile.write(inlines[i])
                    else:
                        outfile.write(inlines[i])
                outfile.close()
    
    def performGroundStateSearch(self, iteration):
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo','\nPerforming ground state search for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                subprocess.call([self.uncleExec, '21'], stdout=self.uncleOut)
                os.chdir(lastDir)

    def makePlots(self, iteration):
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo', '\nMaking plot of ground states and VASP binding energy for ' + atom + '. . .\n'])
                os.chdir(gssDir)                                   
                subprocess.call(['mv', '../vaspFE.out', '.'])
                subprocess.call(['mv', '../vaspBE.out', '.'])
                subprocess.call(['mv', '../vaspHFE.out', '.'])  
                subprocess.call(['gnuplot', 'gss_plot.gp'])
                subprocess.call(['gnuplot', 'BE_plot.gp'])
                subprocess.call(['gnuplot', 'HFE_plot.gp'])
                subprocess.call(['mv','gss.out','gss_' + str(iteration) + '.out'])
                subprocess.call(['mv','gss.pdf','gss_' + str(iteration) + '.pdf'])
        os.chdir(lastDir)
    
    def contains(self, struct, alist):
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
        return False
    
    def getAllGSSStructures(self, iteration, failedStructs):
        allStructs = []
        for n in xrange(len(self.atoms)):
            atomStructs = [] 
            structsByEnergy = []
            gssFile = os.getcwd() + '/' + self.atoms[n] + '/gss/gss_' + str(iteration) + '.out'
            infile = open(gssFile, 'r')
            
            for i,line in enumerate(infile):
                if i >= 2:
                    formEnergy = float(line.strip().split()[7])
                    struct = int(line.strip().split()[0])
                    conc = float(line.strip().split()[1]) #metal concentration
                    # Do not include structures that failed VASP calculations.
                    if not self.contains(struct, failedStructs[n]):
                        structsByEnergy.append([formEnergy, struct])
            infile.close()
            
            structsByEnergy.sort()
            
            for struct in structsByEnergy:
                atomStructs.append(str(struct[1]))
            
            allStructs.append(atomStructs)

                
            
        return allStructs
    
    def getPriorities(self,iteration):  #bch
        lastDir = os.getcwd()
        dir = lastDir
        self.getNcs()
        Ntot = sum(self.Ncs) 
        numberCs = len(self.Ncs)               
        self.priorities = zeros((len(self.atoms),Ntot),dtype = [('struct', 'S10'), ('prior', float)])
#        e_cutoff = zeros(numberCs,dtype = float)
        for iatom in range(len(self.atoms)):
            atomDir = dir + '/' + self.atoms[iatom] + '/gss'
            gss_info = zeros((Ntot),dtype = [('struct', 'S10'), ('conc', float), ('energy', float)]) #columns:  struct, concentration, energy
            lines = self.readfile(atomDir + '/gss_' + str(iteration) + '.out')
            for i,line in enumerate(lines[2:]):
                struct = line.strip().split()[0]
                conc = float(line.strip().split()[2]) #metal concentration
                formEnergy = float(line.strip().split()[7])                
                gss_info[i-2]['struct'] = struct
                gss_info[i-2]['conc'] = conc
                gss_info[i-2]['energy'] = formEnergy
            gss_info = sort(gss_info, order=['conc','energy']) #sorts low to high
            for i in range(Ntot):
                print i, gss_info[i]
            cutoff_percent = 0.05 #approx percentile of energy sorted structs for each concentration
#            de = 0.1 #eV: width of cutoff
            iplace = 0
            for ic,Nc in enumerate(self.Ncs):
#                icutoff = iplace + floor(cutoff_percent*Nc) 
                icutoff = iplace #lowest energy for this concentration
#                ecut = gss_info[icutoff]['energy']
#                print ic, floor(cutoff_percent*Nc), ecut
                width = ceil(cutoff_percent*Nc) 
                for n in range(Nc):
                    istr = iplace+n
                    self.priorities[iatom,istr]['struct'] = gss_info[istr]['struct']

                    en = gss_info[istr]['energy']
                    self.priorities[iatom,istr]['prior'] = 100 * exp(-(istr-icutoff)/width) #100/(exp((istr-icutoff)/width)+1)
                iplace += Nc
            priorfile = open(atomDir+'/priorities.out','w')
            priorfile.write('structure,priority,concentration,FEnergy\n')
            for i in range(Ntot):
                priorfile.write('{:10s}{:10.6f}{:8.4f}{:10.6f}\n'.format( \
                    self.priorities[iatom,i]['struct'],self.priorities[iatom,i]['prior'], \
                       gss_info[i]['conc'],gss_info[i]['energy']))
            priorfile.close()
#            print sort(self.priorities[iatom,:], order = ['prior'] )
                       
            
    def getNcs(self): #bch
        '''Find the number of structures at each concentration''' 
        lastDir = os.getcwd()
        dir = os.getcwd() + '/' + self.atoms[0] + '/gss/'
        os.chdir(dir)
        os.system('sort -g -k 3 gss_1.out > tempout') #sorts by column 3, metal concentration
        lines = self.readfile('tempout');os.system('rm tempout')         
        conc_old = 1.0 #starts with highest concentration first 
        Nc = 0  
        for line in lines[2:]:
            conc = float(line.strip().split()[1])
            if conc == conc_old:
                Nc += 1
            else: 
                self.Ncs.append(Nc) #for the previous concentration           
                Nc = 1
                conc_old = conc
        self.Ncs.append(1) #for the pure H conentration
        os.chdir(lastDir)
        
 
#    def buildPriorityNcs(self,Ncmin,Ntot): #bch
#        '''Ncp in the priority list differs vs concentration. First Ncp=Ncmin, and then Nc is adjusted 
#         at each concentration so those added are roughly proportional to the total number of states 
#         at that concentration (Nc), until all the Ncp's sum to Ntot.  Only need to do this in iteration 1.
#         and for the first atom. Independent of atom type'''  
 
#         
#        gssfile = open('gss.out','r')
            
    def getPriorityStructs(self,Ncmin,Ntot):   #bch
        '''Builds a list of Ntot structures that represent the lowest energy Ncp structures 
        at each concentration'''  
         
         #get number of 
            







        