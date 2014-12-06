import os, subprocess

class RunVasp:
    """ This class is responsible for preparing the directories, retrieving the needed files, and
        submitting VASP jobs to the supercomputer.  It keeps track of the SLURM job ids of all the
        jobs that are currently running from a particular instance of the class. """
    
    def __init__(self, atoms):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.neededFilesFolder = os.getcwd() + '/needed_files/'
        
        self.currJobIds = []

    def addStructName(self,nameadd):
        print os.getcwd()
        jobfile = open('job','r')
        lines = jobfile.readlines()
        print lines
        jobfile.close()
        jobfile = open('job','w')
        for line in lines:
            if 'job-name' in line: 
                jobfile.write(line.strip('\n') + '_'+ nameadd + '\n')
            else:
                jobfile.write(line)
        jobfile.close()  
        
    def clearCurrentJobIds(self):
        self.currJobIds = []
                
    def convergeCheck(self, folder, NSW):
        """Tests whether force convergence is done by whether the last line of Oszicar is less than NSW."""
        try:
            value = self.getSteps(folder)
            return value < NSW #True/False
        except:
            return False #True/False

    def fillDirectories(self, structList):
        """ Fills all the directories with the needed files for VASP to run, namely POSCAR, POTCAR, KPOINTS, 
            INCAR, job, and the VASP executable file. """
            
        for i in xrange(len(self.atoms)):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + self.atoms[i]
            
            os.chdir(atomDir)
            structures = []
            for item in structList[i]:
                if os.path.isdir(item):
                    structures.append(item)
            
            for structure in structures:
                structureDir = os.path.abspath(structure)
                subprocess.call(['cp', 'KPOINTS', 'INCAR', 'job', 'vasp533', structureDir]) # bch
                poscar = open(structureDir + '/POSCAR','r')
                poscarLines = [line.strip() for line in poscar]
                poscar.close()
                
                if poscarLines[0].split()[1] == 'H':
                    subprocess.call(['cp','CH_POTCAR',structureDir + '/POTCAR'])
                elif poscarLines[0].split()[1] == 'M':
                    subprocess.call(['cp','C' + self.atoms[i] + '_POTCAR',structureDir + '/POTCAR'])
                else:
                    subprocess.call(['cp','POTCAR',structureDir])
            
            os.chdir(lastDir)

    def finishCheck(self, folder):
        """ Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR.  The input
        parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """
        lastfolder = os.getcwd()
        os.chdir(os.path.abspath(folder))
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)    
        return newstring[0].find('Voluntary') > -1 #True/False
  
    def getCurrentJobIds(self):
        return self.currJobIds
        
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

    def hasFinished(self, folder):
        if self.finishCheck(folder) and self.convergeCheck(folder, 400):
            return True
        else:
            return False
        
    def makeDOSDirectories(self, structList):
        """ Creates the directories with the needed files for a Density of States run in VASP.
            The directories are created as sub-directories of the original structure directories. """
            
        topDir = os.getcwd()
        for i in xrange(len(self.atoms)):
            elementDir = topDir + '/' + self.atoms[i]
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in structList[i]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        normalDir = structDir + '/normal'
                        if os.path.isdir(normalDir) and self.finishCheck(normalDir) and self.convergeCheck(normalDir, 400):
                            os.chdir(structDir)
                            subprocess.call(['mkdir', 'DOS'])
                            subprocess.call(['cp','normal/CONTCAR','normal/DOSCAR','normal/EIGENVAL',
                                             'normal/IBZKPT','normal/KPOINTS','normal/vasp533',
                                             'normal/OSZICAR','normal/OUTCAR','normal/PCDAT',
                                             'normal/POSCAR','normal/POTCAR','normal/REPORT',
                                             'normal/vasprun.xml','normal/XDATCAR','DOS'])
                            subprocess.call(['mv','normal/CHG','normal/CHGCAR','normal/WAVECAR','DOS'])#bch                     

                            self.makeDOS_INCAR()
                            self.makeDOSJobFile(self.atoms[i]+structure) #bch 
                            subprocess.call(['cp','DOS/CONTCAR','DOS/POSCAR'])
                            os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir)
    
    def makeDOS_INCAR(self):
        """ Creates an INCAR file for the Density of States run in VASP.  The notable changes are:
                IBRION=-1 -- This tells VASP not to move the ions.
                NSW=0     -- This tells VASP that there will be no ionic relaxation steps.
                LORBIT=11 -- This creates the PROCAR file which can be used to project onto the C, H, and M atoms. """
                
        incar = open('DOS/INCAR','w')
        
        incar.write("IBRION=-1\n")
        incar.write("ISIF=4\n")
        incar.write("NSW=0\n")
        incar.write("Algo=Normal\n")
        incar.write("PREC=Normal\n")
        incar.write("EDIFF=5E-4\n")
        incar.write("EDIFFG=5E-4\n")
        incar.write("ISMEAR=0\n")
        incar.write("ISPIN=2\n")
        incar.write("LREAL=Auto\n")
        incar.write("LORBIT=11\n")
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")
        incar.write("NPAR = 4\n")    
        incar.close()

    def makeDOSJobFile(self,name):
        
        jobFile = open('DOS/job','w')
        
        jobFile.write("#!/bin/bash\n\n")
        jobFile.write("#SBATCH --time=6:00:00\n")
        jobFile.write("#SBATCH --ntasks=16\n")
        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
        jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
        jobFile.write("#SBATCH --mail-type=FAIL\n")
        jobFile.write("#SBATCH --mail-type=END\n")  
        jobFile.write("#SBATCH --job-name=%s\n" % name)
        jobFile.write("\nmpiexec vasp533 > vasp.out\n")    
        jobFile.close() 
            
    def makeJobFiles(self):#bch name
        """ Creates a standard job file for submitting a VASP job on the supercomputer. """
    
        dirList = self.atoms
        
        for direc in dirList:
            name = direc #bch (which atom)
            jobFile = open(direc + '/job','w')
    
            jobFile.write("#!/bin/bash\n\n")
            jobFile.write("#SBATCH --time=06:00:00\n")
            jobFile.write("#SBATCH --ntasks=16\n")
            jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
            jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")              
            jobFile.write("#SBATCH --mail-type=FAIL\n")
            jobFile.write("#SBATCH --mail-type=END\n") 
            jobFile.write("#SBATCH --job-name=%s\n" % name) #bch  adds atom name.  Later we add structure name                    
            jobFile.write("\nmpiexec vasp533 > vasp.out\n")
    
            jobFile.close()    
    
    def makeKPOINTS(self, num1, num2):
        """ Creates a KPOINTS file based on the input parameters num1 and num2. It specifies that the job will have
        num1 x num2 kpoints. For example, if we wanted to specify an 8x8 kpoints scenario, we would call
        makeKPOINTS(8, 8). """
        
        dirList = self.atoms
        
        for direc in dirList:
            kpoints = open(direc + '/KPOINTS','w')
    
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write(str(num1) + ' ' + str(num2) + ' 1\n')
            kpoints.write('0 0 0')
    
            kpoints.close()
   
    def makeLowINCARs(self):
        """ Creates a standard INCAR file and puts it in each different structure's top directory. """
        
        dirList = self.atoms
        
        for direc in dirList:
            incar = open(direc + '/INCAR','w')
    
            incar.write("IBRION=2\n")
            incar.write("ISIF=4\n")
            incar.write("NSW=400\n")
            incar.write("Algo=VeryFast\n")
            incar.write("PREC=Low\n") #BCH was LOW!
            incar.write("EDIFF=2E-3\n")
            incar.write("EDIFFG=2E-3\n")
            incar.write("ISMEAR=0\n")
            incar.write("ISPIN=2\n")
            incar.write("LREAL=Auto\n")
            incar.write("SIGMA=0.1\n")
            incar.write("LWAVE=.TRUE.\n")
            incar.write("LCHARG=.TRUE.\n")
            incar.write("NPAR = 4\n") 
            #incar.write("BMIX = 0.3\n") #BCH 
            incar.close()
    
    def makeNormalDirectories(self, structList):
        topDir = os.getcwd()
        for i in xrange(len(self.atoms)):
            elementDir = topDir + '/' + self.atoms[i]
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in structList[i]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir) and self.finishCheck(structDir) and self.convergeCheck(structDir, 400):
                        os.chdir(structDir)
                        subprocess.call(['mkdir', 'normal'])
                        subprocess.call(['cp','CONTCAR','DOSCAR','EIGENVAL',
                                         'IBZKPT','KPOINTS','vasp533',
                                         'OSZICAR','OUTCAR','PCDAT',
                                         'POSCAR','POTCAR','REPORT',
                                         'vasprun.xml','job','XDATCAR','normal'])
                        subprocess.call(['mv','CHG','CHGCAR','WAVECAR','normal'])#bch                     
                        self.makeNormalINCAR()
                        subprocess.call(['cp','normal/CONTCAR','normal/POSCAR'])
                        os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])            
            os.chdir(topDir)

    def makeNormalINCAR(self):
        """ Creates a standard INCAR file and puts it in each different structure's top directory. """
        
        incar = open('normal/INCAR','w')
    
        incar.write("IBRION=2\n")
        incar.write("ISIF=4\n")
        incar.write("NSW=400\n")
        incar.write("Algo=VeryFast\n")
        incar.write("PREC=Normal\n")
        incar.write("EDIFF=5E-4\n")
        incar.write("EDIFFG=5E-4\n")
        incar.write("ISMEAR=0\n")
        incar.write("ISPIN=2\n")
        incar.write("LREAL=Auto\n")
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")
        incar.write("NPAR = 4\n")   
        #incar.write("BMIX = 0.3\n") #BCH          
        incar.close()

    def makePOTCARs(self):
        """ Creates a POTCAR file based on the atoms in the given input list. Concatenates the individual POTCAR
        files to make a single POTCAR file for the multi-atom structure. It is assumed that the atoms in the
        list are given in the correct order. """
    
        for atom in self.atoms:
            CPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/C/POTCAR"
            HPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/H/POTCAR"
            atomPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + "/POTCAR"
            if os.path.exists(atomPotcarDir):
                CPotcar = open(CPotcarDir, 'r')
                CLines = CPotcar.readlines()
                CPotcar.close()
                    
                HPotcar = open(HPotcarDir, 'r')
                HLines = HPotcar.readlines()
                HPotcar.close()
                    
                atomPotcar = open(atomPotcarDir,'r')
                atomLines = atomPotcar.readlines()
                atomPotcar.close()
            
                potcar = open(atom + '/POTCAR', 'w')
                for line in CLines:
                    potcar.write(line)
                        
                for line in HLines:
                    potcar.write(line)
                    
                for line in atomLines:
                    potcar.write(line)
                        
                potcar.close()
                    
            else:
                subprocess.call(['echo','ERROR: Could not find a POTCAR file for \'' + atom + '\''])
                subprocess.call(['echo','Removing POTCAR . . .'])
                potcar.close()
                subprocess.call(['rm','POTCAR'])
                return # ? indent?  bch
    def makePurePOTCARs(self):
        """ Some of the structures that need to be submitted to VASP for relaxation are what we call
            "pure" structures.  This means that (other than carbon atoms) the structure only contains 
            one other type of atom.  In the binary representation of the structure, this means that it
            is either all '1's or all '0's.  VASP gives a segmentation fault even if we tell it that 
            there are zero of one of the kinds of atoms.  It needs a POTCAR file that doesn't even mention
            the element that is not a part of the "pure" structure.  This method creates these POTCAR 
            files. """
        
        for atom in self.atoms:
            atomPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom
            
            if os.path.isdir(atomPotcarDir):
                purePotcar = open(atom + "/C" + atom + "_POTCAR", "w")
                
                CPotcar = open("/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/C/POTCAR", "r")
                CLines = CPotcar.readlines()
                CPotcar.close()
                
                atomPotcar = open(atomPotcarDir + "/POTCAR", "r")
                atomLines = atomPotcar.readlines()
                atomPotcar.close()
                
                for line in CLines:
                    purePotcar.write(line)
                
                for line in atomLines:
                    purePotcar.write(line)
                
                purePotcar.close()
                
                CHPotcar = open(atom + "/CH_POTCAR", "w")
                
                HPotcar = open("/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/H/POTCAR", "r")
                HLines = HPotcar.readlines()
                HPotcar.close()
                
                for line in CLines:
                    CHPotcar.write(line)
                
                for line in HLines:
                    CHPotcar.write(line)
                
                CHPotcar.close()    
                
    def makeRunHexMono(self): #bch all
        topDir = os.getcwd()
        if not os.path.isdir('hex_monolayer_refs'): os.mkdir('hex_monolayer_refs')
        os.chdir('hex_monolayer_refs')
#        os.system('rm -r -f *')
        for atom in self.atoms:
            os.mkdir(atom)
            os.chdir(atom)
            atomPotcar = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + '/POTCAR'
            if os.path.exists(atomPotcar):
                subprocess.call(['cp', atomPotcar, '.'])
            else:
                system.exit('Failed to read POTCAR in makeSingleDirectories')                
            jobFile = open('job','w')
            jobFile.write("#!/bin/bash\n\n")
            jobFile.write("#SBATCH --time=03:00:00\n")
            jobFile.write("#SBATCH --ntasks=1\n")
            jobFile.write("#SBATCH --mem-per-cpu=4G\n")
            jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")              
            jobFile.write("#SBATCH --mail-type=FAIL\n")
            jobFile.write("#SBATCH --mail-type=END\n") 
            jobFile.write("#SBATCH --job-name=hexm%s\n" % atom)           
            jobFile.write("\nmpiexec vasp533 > vasp.out\n") 
            jobFile.close()
            incar = open('INCAR','w')
            incar.write("IBRION=2\n")
            incar.write("ISIF=4\n")
            incar.write("NSW=400\n")
            incar.write("PREC=High\n")
            incar.write("EDIFF=1E-6\n")
            incar.write("ISPIN=2\n")
            incar.write("LWAVE=.FALSE.\n")
            incar.write("LCHARG=.FALSE.\n")         
            incar.close()            
            kpoints = open('KPOINTS','w')
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write('8 8 8\n')
            kpoints.write('0 0 0')
            kpoints.close()
            poscar = open('POSCAR','w') 
            poscar.write('Hexagonal metal monolayer\n')
            poscar.write('1.0\n') #scale same as graphene for now
            poscar.write('2.13128850  -1.23050000  0\n') 
            poscar.write('2.13128850   1.23050000 0\n')
            poscar.write('0.00000000   0.00000000  15.00000000\n')
            poscar.write('1\n')
            poscar.write('Cartesian\n')
            poscar.write('0 0 0\n')
            poscar.close()
            subprocess.call(['sbatch','job'])
            os.chdir('../')      
        os.chdir(topDir)   
        
    def makeRunSingleDirectories(self): #bch all
        topDir = os.getcwd()
        if not os.path.isdir('single_atoms'): os.mkdir('single_atoms')
        os.chdir('single_atoms')
#        os.system('rm -r -f *')
        for atom in self.atoms:
            os.mkdir(atom)
            os.chdir(atom)
            atomPotcar = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + '/POTCAR'
            if os.path.exists(atomPotcar):
                subprocess.call(['cp', atomPotcar, '.'])
            else:
                system.exit('Failed to read POTCAR in makeSingleDirectories')                
            jobFile = open('job','w')
            jobFile.write("#!/bin/bash\n\n")
            jobFile.write("#SBATCH --time=03:00:00\n")
            jobFile.write("#SBATCH --ntasks=1\n")
            jobFile.write("#SBATCH --mem-per-cpu=4G\n")
            jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")              
            jobFile.write("#SBATCH --mail-type=FAIL\n")
            jobFile.write("#SBATCH --mail-type=END\n") 
            jobFile.write("#SBATCH --job-name=isol_%s\n" % atom)           
            jobFile.write("\nmpiexec vasp533 > vasp.out\n") 
            jobFile.close()
            incar = open('INCAR','w')
            incar.write("IBRION=-1\n")
            incar.write("NELM=400\n")
            incar.write("PREC=High\n")
            incar.write("EDIFF=1E-6\n")
            incar.write("ISPIN=2\n")
            incar.write("LWAVE=.FALSE.\n")
            incar.write("LCHARG=.FALSE.\n")         
            incar.close()            
            kpoints = open('KPOINTS','w')
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write('1 1 1\n')
            kpoints.write('0 0 0')
            kpoints.close()
            poscar = open('POSCAR','w') 
            poscar.write('Isolated atom\n')
            poscar.write('1\n')
            poscar.write('20 0 0\n') 
            poscar.write('0 20 0\n')
            poscar.write('0.0 0.0 20.0\n')
            poscar.write('1\n')
            poscar.write('Cartesian\n')
            poscar.write('0 0 0\n')
            poscar.close()
            subprocess.call(['sbatch','job'])
            os.chdir('../')      
        os.chdir(topDir)  
                     
    def prepareForVasp(self, structList):
        self.makeLowINCARs()
        self.makePurePOTCARs()
        self.makePOTCARs()
        self.makeKPOINTS(6, 6) #bch changed from 8 8
        self.makeJobFiles() 
        self.linkVaspExec()       
        self.fillDirectories(structList)
                                                  
    def run(self, runNum, structList):       
        if runNum == 1:
            self.startJobs(structList)
    
        elif runNum == 2:
            self.makeNormalDirectories(structList)
            self.startNormalJobs(structList)
           
        elif runNum == 3:
            self.makeDOSDirectories(structList)
            self.startDOSJobs(structList)

    def startDOSJobs(self, structList):
        """ Submits all the VASP jobs to the supercomputer for the Density of States run. """
        
        topDir = os.getcwd()
        self.clearCurrentJobIds()
        
        for i in xrange(len(self.atoms)):
            elementDir = topDir + '/' + self.atoms[i]
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                
                for structure in structList[i]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        os.chdir(structDir)
                        
                        dosDir = structDir + '/DOS'
                        if os.path.isdir(dosDir):
                            os.chdir(dosDir)
                            proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                            jobid = proc.communicate()[0].split()[3]
                            subprocess.call(['echo','Submitted job ' + jobid])
                            self.currJobIds.append(jobid)
                        
                        os.chdir(structDir)
                    
                    os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir)
            
    def startJobs(self, structList):
        """ Submits all the VASP jobs to the supercomputer for initial ionic relaxation. """
        
        self.clearCurrentJobIds()
        
        for i in xrange(len(self.atoms)):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + self.atoms[i]
            
            os.chdir(atomDir)
            
            structures = []
            for item in structList[i]:
                if os.path.isdir(item):
                    structures.append(item)
            
            for structure in structures:
#                if not(self.finishCheck(structure) and self.convergeCheck(structure, 2)): #finalDir
                    #bch put this case here to not run finished jobs
                os.chdir(structure)
                self.addStructName(structure) #bch adds structure to name
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
                subprocess.call(['echo', 'Submitted job ' + jobid])
                self.currJobIds.append(jobid)
                os.chdir(atomDir)
            
            os.chdir(lastDir)
                  
    def startNormalJobs(self, structList):
        self.clearCurrentJobIds()
        
        for i in xrange(len(self.atoms)):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + self.atoms[i]
            
            os.chdir(atomDir)
            
            structures = []
            for item in structList[i]:
                if os.path.isdir(item + '/normal'):
                    structures.append(item)
            
            for structure in structures:
                os.chdir(structure + '/normal')
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
                subprocess.call(['echo','Submitted job ' + jobid])
                self.currJobIds.append(jobid)
                os.chdir(atomDir)
            
            os.chdir(lastDir)                
        