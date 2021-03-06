'''
Created on Aug 29, 2014


'''
from numpy import zeros, mod #bch
import os, subprocess
from random import random


class MakeUncleFiles:


    def __init__(self, atoms):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.structuresInLengths = zeros(len(self.atoms))
    
        self.structList = []
        self.newlyFailed = []
        self.pureHenergy = 0.0
        self.pureMenergy = 0.0
        
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0
        self.singleE = [] #bch
        self.hexE = [] #bch
        self.finalDir = '' #bch
    

    def closeOutFiles(self):
        """ Closes both the structures.in and structures.holdout files. """
        self.infile.close()
#        self.holdoutFile.close()

    def convergeCheck(self, folder, NSW):
        """Tests whether ionic convergence is done, by whether the last line of Oszicar is less than NSW."""
#        try:
        value = self.getSteps(folder)
        return value < NSW and self.energyDropCheck(folder)
#        except:
#            return False  

    def elConvergeCheck(self,folder,NELM): #bch
        """Tests electronic convergence is done by whether the electronic step is less than NELM."""
        try:
            value = self.getElSteps(folder)
            return value < NELM #True/False
        except:
            return False #True/False
        
    def electronicConvergeFinish(self,dir): #bch
        '''Test requires electronic convergence AND vasp finishing'''
        #get NELM, the max electronic steps allowed in the run.  Using first directory in dirslist
        proc = subprocess.Popen(['grep','-i','NELM',dir+'/INCAR'],stdout=subprocess.PIPE)
        result =  proc.communicate()[0]
        NELM = int(result.split('=')[1].split()[0])
        return self.elConvergeCheck(dir,NELM) and finishCheck(dir)
    
    def energyDropCheck(self,dir):
        '''tests whether the energies have dropped in OSZICAR...rising energies are unphysical 
        and show a numerical convergence problem. The factor like 0.99 allows for very small energy rises only'''
        lines = readfile(dir+'/OSZICAR') 
        energies = []
        for line in lines:
            if 'F=' in line:
                energies.append(float(line.split()[2]))
        return energies[-1] <= 0.99*energies[0]    
   
            
        lastfolder = os.getcwd()
        os.chdir(dir) 
        proc = subprocess.Popen(['grep', 'Voluntary', dir+'/OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)   
        return newstring[0].find('Voluntary') > -1 #True/False
    

    def getElSteps(self,folder):#bch
        '''number of electronic steps for isolated runs'''
        lastfolder = os.getcwd()
        os.chdir(folder)
        try:
            oszicar = open('OSZICAR','r') 
            laststep = oszicar.readlines()[-2].split()[1] # Next to last line, 2nd entry
            oszicar.close()
            os.chdir(lastfolder) 
            value = int(laststep)
            return value         
        except:
            os.chdir(lastfolder)         
            return 9999    
        
    def getEnergy(self,dir): #bch
        lines = readfile(dir+'/OSZICAR')
        if len(lines[-1].split())>1:
            energy = float(lines[-1].split()[2])  #Free energy
        else: 
            energy = 0.0
        return energy

    def getNSW(self,dir): #bch
        proc = subprocess.Popen(['grep','-i','NSW',dir+'/INCAR'],stdout=subprocess.PIPE) 
        return int(proc.communicate()[0].split('=')[-1])   

    def getSteps(self, folder):
        '''number of steps in relaxation, as an integer'''
        lastfolder = os.getcwd()
        os.chdir(folder)
        if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
            os.chdir(lastfolder) 
            return -9999
        oszicar = open('OSZICAR','r')
        laststep = oszicar.readlines()[-1].split()[0]
        oszicar.close()
        os.chdir(lastfolder)  
        try:
            value = int(laststep)
            return value
        except:
            return 9999   

    def getStructuresInLengths(self):
        lengths = zeros(len(self.structuresInLengths))
        for i in xrange(len(self.structuresInLengths)):
            lengths[i] = self.structuresInLengths[i]
        
        return lengths

    def getStructureList(self):
        """ Returns the list of usable structures. """
        returnList = []
        for i in xrange(len(self.structList)):
            subList = []
            for j in xrange(len(self.structList[i])):
                subList.append(self.structList[i][j].split('/')[-1])
            returnList.append(subList)
        
        failedList = []
        for i in xrange(len(self.newlyFailed)):
            subList = []
            for j in xrange(len(self.newlyFailed[i])):
                subList.append(self.newlyFailed[i][j].split('/')[-1])
            failedList.append(subList)
        
        return returnList, failedList
           
    def hexMonolayerEnergies(self,dir1): #bch
        file = open(dir1 +'/hex_monolayer_refs/hex_energies','w')
        self.hexE = zeros(len(self.atoms),dtype = float) +100  #default to large number so can tell if not read
        subprocess.call(['echo', '\nReading hexagonal monolayer energies\n'])
        for i,atom in enumerate(self.atoms):
            dir2 = dir1 + '/hex_monolayer_refs'+'/'+atom
            if finishCheck(dir2) and convergeCheck(dir2, self.getNSW(dir2)): #finalDir
                print'{} monolayer (per atom): {:8.4f} '.format(atom,self.getEnergy(dir2))
                file.write('{} monolayer (per atom): {:8.4f} \n'.format(atom,self.getEnergy(dir2)))
                self.hexE[i] = self.getEnergy(dir2) 
            else:
                file.write('{} monolayer not converged \n'.format(atom))
        os.chdir(dir1)               

    def makeUncleFiles(self):
        self.singleAtomsEnergies(os.getcwd())  #bch
        self.hexMonolayerEnergies(os.getcwd())   #bch
        
        self.setStructureList()
       
        for i in xrange(len(self.atoms)):
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                subprocess.call(['echo', '\nCreating structures.in and structures.holdout files for ' + self.atoms[i] + '\n'])
                self.initialize()
                self.infile = open(atomDir + '/structures.in','w')
                self.holdoutFile = open(atomDir + '/structures.holdout','w')
                self.sortStructsByFormEnergy(i)
                self.writeHeader()        
                num = 0
                structuresInCount = 0
                for structure in self.structList[i]:
#                    if structuresInCount >= 500:    # Write a maximum of 500 structures to the file
#                        break                       # for any given atom.
                    whichFile = self.writeUnclePOSCAR(structure, i)
                    if whichFile == 'in':
                        structuresInCount += 1
                    num += 1
                
                self.structuresInLengths[i] = structuresInCount
                self.closeOutFiles()
             
    def readfile(self,filepath): #bch
        file1 = open(filepath,'r')
        lines = file1.readlines()
        file1.close()
        return lines
    
    def initialize(self):
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0  


    def setAtomCounts(self, poscarDir):
        """ Retrieves the number of H atoms and the number of M atoms from the POSCAR file
            and sets the corresponding members. """
        self.atomCounts = []

        poscar = open(poscarDir + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        counts = poscarLines[5].strip().split()
        
        if len(counts) == 3:
            self.atomCounts.append(int(counts[1]))
            self.atomCounts.append(int(counts[2]))
        elif len(counts) == 2:
            if poscarLines[0].split()[1] == 'H':
                self.atomCounts.append(int(counts[1]))
                self.atomCounts.append(0)
            elif poscarLines[0].split()[1] == 'M':
                self.atomCounts.append(0)
                self.atomCounts.append(int(counts[1]))

    def setAtomPositions(self, poscarDir):
        """ Retrieves the positions of each of the atoms.  Appends the x-coordinate to the xPositions
            list, the y-coordinate to the yPositions list.  For a surface in UNCLE the z-coordinate
            is always zero. """
        poscar = open(poscarDir + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomPositions = poscarLines[7 + sum(self.atomCounts):7 + (2 * sum(self.atomCounts))]
        self.atomPositions = [line.strip().split() for line in self.atomPositions]
           
        for pos in self.atomPositions:
            self.xPositions.append(float(pos[0]))
            self.yPositions.append(float(pos[1]))
            self.zPositions.append(0.0)

    def setEnergy(self, directory):  
        """ Retrieves the energy of the structure from the OSZICAR file and sets the corresponding 
            member. """   
        try:
            oszicar = open(directory + self.finalDir + '/OSZICAR','r')
            energy = oszicar.readlines()[-1].split()[2]
            oszicar.close()
        except:
            energy = 0
        
        energy = float(energy)
        peratom = energy / sum(self.atomCounts)
        
        self.energy = str(peratom)

    def setIDString(self, poscarDir):
        """ Sets the first written line of each structure to the form:
                C H (Metal Atom)  Structure:  #(Decimal) (#(Binary)) """
        poscar = open(poscarDir + '/POSCAR', 'r')
        ID = poscar.readlines()[0].strip()
        
        self.idString = ID

    def setLatticeVectors(self, structFile):
        """ Gets the lattice vectors from the first structure in the structList and sets
            the corresponding member components. """
        vecFile = open(structFile + '/POSCAR','r')
        vecFileLines = vecFile.readlines()
        vecFile.close()
        
        vec1 = vecFileLines[2].strip().split()
        vec2 = vecFileLines[3].strip().split()
        vec3 = vecFileLines[4].strip().split()
        
        vec1comps = [float(comp) for comp in vec1]
        vec2comps = [float(comp) for comp in vec2]
        vec3comps = [float(comp) for comp in vec3]
            
        
        self.vec1x = vec1comps[0]
        self.vec1y = vec1comps[1]
        if vec1comps[2] == 15.0:
            self.vec1z = 1000.0
        else:
            self.vec1z = vec1comps[2]
        
        self.vec2x = vec2comps[0]
        self.vec2y = vec2comps[1]
        if vec2comps[2] == 15.0:
            self.vec2z = 1000.0
        else:
            self.vec2z = vec2comps[2]
        
        self.vec3x = vec3comps[0]
        self.vec3y = vec3comps[1]
        if vec3comps[2] == 15.0:
            self.vec3z = 1000.0
        else:
            self.vec3z = vec3comps[2]
    
    def setStructureList(self):
        """ Initializes the list of structures to add to the structures.in and structures.holdout
            files by adding only the structures that VASP finished relaxing to the member
            structList. Sorts the list by concentration. """
        
        self.structList = []
        
        lastDir = os.getcwd()
        for atom in self.atoms:
            atomDir = lastDir + '/' + atom
            os.chdir(atomDir)
            pureHdir = os.getcwd() + '/1'
            pureMdir = os.getcwd() + '/3'
        
            self.setAtomCounts(pureHdir)
            self.setEnergy(pureHdir)
            self.pureHenergy = float(self.energy)
        
            self.setAtomCounts(pureMdir)
            self.setEnergy(pureMdir)
            self.pureMenergy = float(self.energy)
        
            conclist = []
            atomStructs = []
            failed = []
            dirList = os.listdir(atomDir)            
            dirList2 = [] #bch
            for i in dirList:#bch
                if os.path.isdir(i): dirList2.append(i) #bch
            for i, item in enumerate(dirList2):#bch enumerate
                if mod(i+1,100) == 0: print 'Checking',i+1,'of',len(dirList2), 'structures in', atom  #bch
                fullPath = os.path.abspath(item)
                if os.path.isdir(fullPath):
#                if os.path.isdir(self.finalDir): #bch finalDir indtead of '/DOS'
#                        print finishCheck(fullPath + self.finalDir) , convergeCheck(fullPath + self.finalDir, self.getNSW()) #BCH nsw
                    if finishCheck(fullPath + self.finalDir) and convergeCheck(fullPath + self.finalDir, self.getNSW(fullPath + self.finalDir)): #finalDir
                        # Check for concentration
                        self.setAtomCounts(fullPath)
                        concentration = 0.0
                        if self.atomCounts[0] == 0:
                            concentration = 1.0
                        else:
                            concentration = float(float(self.atomCounts[1]) / float(self.atomCounts[0] + self.atomCounts[1]))
                    
                        conclist.append([concentration, fullPath])
                    else:
                        failed.append(fullPath)
                else:
                    # Don't add the 'gss' or 'fits' directories.
                    if not fullPath.split('/')[-1] == 'gss' and not fullPath.split('/')[-1] == 'fits':
                        failed.append(fullPath)
    
            self.newlyFailed.append(failed)
            conclist.sort()
            
            for i in xrange(len(conclist)):
                atomStructs.append(conclist[i][1])
            self.structList.append(atomStructs)
                        
            os.chdir(lastDir)

    def singleAtomsEnergies(self,dir1): #bch
        self.singleE = zeros(len(self.atoms),dtype = float) +100  #default to large number so can tell if not read
        subprocess.call(['echo', '\nReading single atom energies\n'])
        file = open(dir1 +'/single_atoms/single_atom_energies','w')
        for i,atom in enumerate(self.atoms):
            dir2 = dir1 + '/single_atoms'+'/'+atom
            if self.electronicConvergeFinish(dir2): 
                print 'Energy of {} atom: {:8.4f} \n'.format(atom,self.getEnergy(dir2))
                file.write('{} atom: {:12.8f} \n'.format(atom,self.getEnergy(dir2)))
                self.singleE[i] = self.getEnergy(dir2)
        file.close()  
        os.chdir(dir1)   
            
    def sortStructsByFormEnergy(self, atomInd):
        '''Besides sorting, writes vasp formation energies vs concentration to vaspFE.out so can be plotted with fitted data.
        Also writes vaspBE, which gives the structure's binding energy per adatom relative to graphene and separated adatoms
        Also writes vaspHFE, which is a formation energy with the metal atom taken from a relaxed hexagonal monolayer of the metal'''  #bch  
             
        lastDir = os.getcwd()
        os.chdir(lastDir + '/' + self.atoms[atomInd])
        pureHdir = os.getcwd() + '/1'
        pureMdir = os.getcwd() + '/3'
        
        self.setAtomCounts(pureHdir)
        self.setEnergy(pureHdir)
        self.pureHenergy = float(self.energy)
        
        self.setAtomCounts(pureMdir)
        self.setEnergy(pureMdir)
        self.pureMenergy = float(self.energy)
        eIsolatedH = -1.115 #bch
        eIsolatedC = -1.3179 #bch
        eH2 = -6.7591696/2.0 #bch
        energyGraphene = -18.456521 #for 2 C atoms #bch
        
        formEnergyList = []
        sortedStructs = []
        vaspFEfile = open('vaspFE.out','w') #bch 
        vaspBEfile = open('vaspBE.out','w') #bch 
        vaspHFEfile = open('vaspHFE.out','w') #bch  
        for structDir in self.structList[atomInd]:
            self.setAtomCounts(structDir)
            self.setEnergy(structDir)
            structEnergy = float(self.energy)
            struct = structDir.split('/')[-1] #bch
            concentration = 0.0
            if self.atomCounts[0] == 0:
                concentration = 1.0
            else:
                concentration = float(float(self.atomCounts[1]) / float(self.atomCounts[0] + self.atomCounts[1]))
                        
            formationEnergy = structEnergy - (concentration * self.pureMenergy + (1.0 - concentration) * self.pureHenergy)
            formEnergyList.append([formationEnergy, structDir])
            vaspFEfile.write('{:10s} {:12.8f} {:12.8f}\n'.format(struct,concentration,formationEnergy))#bch 
            
            ncarbon = self.atomCounts[0] + self.atomCounts[1] #bch:  
            bindEnergy = structEnergy - (self.atomCounts[0]*eIsolatedH + self.atomCounts[1]*self.singleE[atomInd] + ncarbon*energyGraphene/2)/ float(self.atomCounts[0] + self.atomCounts[1]) #2 atoms in graphene 
            vaspBEfile.write('{:10s} {:12.8f} {:12.8f}\n'.format(struct,concentration,bindEnergy))#bch  
            hexFormationEnergy = structEnergy - energyGraphene/2  - (concentration * self.hexE[atomInd] + (1.0 - concentration) * eH2)
            vaspHFEfile.write('{:10s} {:12.8f} {:12.8f}\n'.format(struct,concentration,hexFormationEnergy))#bch    
                                      
        vaspFEfile.close()#bch
        vaspBEfile.close()#bch
        vaspHFEfile.close()#bch
        formEnergyList.sort()
        
        for pair in formEnergyList:
            sortedStructs.append(pair[1])
        
        self.structList[atomInd] = sortedStructs
            
        os.chdir(lastDir)
    

    def writeAtomCounts(self):
        """ Writes the number of H atoms and the number of M atoms to the structures.in or 
            structures.holdout file. """
        if len(self.atomCounts) == 2:
            self.outfile.write(str(self.atomCounts[0]) + " " + str(self.atomCounts[1]) + "\n")
        elif len(self.atomCounts) == 1:
            self.outfile.write(str(self.atomCounts[0]) + "\n")
    
    def writeAtomPositions(self):
        """ Writes the positions of the atoms in the current structure to the structures.in
            or structures.holdout file.  The positions are written in the form:
                z-coord   x-coord   y-coord
            because this is the convention that UNCLE uses. """
        self.outfile.write("Cartesian\n")
        for i in xrange(len(self.atomPositions)):
            self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % 
                               (self.zPositions[i], self.xPositions[i], self.yPositions[i]))
    def writeDashedLine(self):
        """ Writes a dashed line in the structures.in/.holdout files as a separator between
            different structures. """
        self.outfile.write("#------------------------------------------------\n")

    def writeEnergy(self):
        """ Writes the energy of the current structure to the structures.in or structures.holdout
            file. """
        self.outfile.write("#Energy:\n")
        self.outfile.write(str(self.energy) + "\n") 

    def writeHeader(self):
        """ Writes the headers of the structures.in and structures.holdout files. """
        self.infile.write(self.header)
#        self.holdoutFile.write(self.header)  

    def writeIDString(self):
        """ Writes the ID string of the current structure to either the structures.in or 
            structures.holdout file. """
        concentration = float(float(self.atomCounts[1]) / float(sum(self.atomCounts)))
        formationEnergy = float(self.energy) - (concentration * self.pureMenergy + (1.0 - concentration) * self.pureHenergy)
        
        self.outfile.write(self.idString + " FE = " + str(formationEnergy) + ", Concentration = " + str(concentration) + "\n")
        
    def writeLatticeVecs(self):
        """ Writes the lattice vectors of the current structure to the structures.in or
            structures.holdout file. """
        self.outfile.write("1.0\n")
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec1z, self.vec1x, self.vec1y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec2z, self.vec2x, self.vec2y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec3z, self.vec3x, self.vec3y))
    
    def writeUnclePOSCAR(self, poscarDir, atomInd):
        """ Calls all the methods needed to write all the needed information about the current
            structure to the structures.in or structures.holdout files.  Puts a maximum of 10%
            of the structures in the structures.holdout file. """
#        if self.holdoutCount / float(len(self.structList[atomInd])) < .10 and random() < .15:
#            self.outfile = self.holdoutFile
#            self.holdoutCount += 1
#        else:
        self.outfile = self.infile #bch removed holdout info
        self.inCount += 1
        
        #self.outfile = self.infile
        self.setIDString(poscarDir)
        self.setLatticeVectors(poscarDir)
        self.setAtomCounts(poscarDir)
        self.setAtomPositions(poscarDir)
        self.setEnergy(poscarDir)
        
        # Make sure the pure structures go in structures.in
        if self.idString.split()[0] == 'PURE':
            self.outfile = self.infile
        
        self.writeDashedLine()
        self.writeIDString()
        self.writeLatticeVecs()
        self.writeAtomCounts()
        self.writeAtomPositions()
        self.writeEnergy()
        
        if self.outfile.name.split('.')[-1] =='holdout':
            return 'holdout'
        else:
            return 'in'
        











