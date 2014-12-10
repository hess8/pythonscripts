import os, subprocess, copy, string  

def nstrip(list):
    '''Strips off /n'''
    import string
    list2 = []
    for string1 in list:   
        string2 = string1.strip("\n")
        list2.append(string2)
    return list2

class VaspTools:

    def __init__(self,mainFolder,runname,runtype,poscar,kpoints,incar,potcar,
        poscarvar,kpointsvar,incarvar,potcarvar,potcardir,toCheckList, 
        checkedList, toRunList ):
        self.MainFolder = mainFolder
        self.RunName = runname
        self.RunType = runtype
        self.Poscar = poscar
        self.PoscarVars = poscarvar
        self.Kpoints = kpoints
        self.KpointsVars = kpointsvar
        self.Incar = incar
        self.IncarVars = incarvar
        self.Potcar = potcar
        self.PotcarVars = potcarvar
        self.PotcarDir = potcardir    
        self.FolderList = []
        self.TotalRuns = 0
        self.ToCheckList = toCheckList
        self.CheckedList = checkedList
        self.ToRunList = toRunList
       
    def AddToList(self,folder):
        files = os.listdir(folder)
        for path in files:
            if os.path.isdir(folder+path+"/"):
                self.ToCheckList.append(folder+path+"/")
                self.AddToList(folder+path+"/")        
    
    def BuildNewRun(self):
        os.chdir(self.MainFolder)
        for folder in self.RunType:
            try:
                os.chdir(folder)
            except OSError:
                os.system('mkdir %s' % folder)
                os.chdir(folder)
        self.DescendPotcar(self.PotcarVars) 
        print "Total Number of run folders created:"
        print self.TotalRuns
        
    def DescendPotcar(self,varmap): #varmap is elementlist. Creates run folder if needed
        varmap = copy.deepcopy(varmap)        
        if not (os.path.isfile('POTCAR')):
            try:
                os.system('cp %s %s' % (self.Potcar,'POTCAR' ))
            except:
                print "File system error trying to create Potcar"               
        if len(varmap) == 0:
            self.DescendPoscar(self.PoscarVars) 
            return           
        key = varmap.items()[0][0]
        label = key[1:] #keyword @adatom        
        varlist = varmap[key]
        del varmap[key]
        for option in varlist:
            try:
                os.chdir(str(label) + '_' + str(option))
            except OSError:
                os.system('mkdir %s' % str(label) + '_' + str(option))
                os.chdir(str(label) + '_' + str(option))
            try:
                os.system("cp ../POTCAR .")
            except:
                print 'File system error'
            self.AlterPotcar(option,key)
            self.DescendPotcar(copy.deepcopy(varmap))
            os.chdir('..')
        try:
            os.system('rm POTCAR')
        except:
            print "File system error\n"
                            
    def DescendPoscar(self,varmap):
        varmap = copy.deepcopy(varmap)
        
        if not (os.path.isfile('POSCAR')):
            try:
                os.system('cp %s %s' % (self.Poscar,'POSCAR' ))
            except:
                print "File system error"                
        if len(varmap) == 0:
            self.DescendKpoints(self.KpointsVars)
            return            
        key = varmap.items()[0][0]
        label = key[1:]
        varlist = varmap[key]
        del varmap[key]        
        for option in varlist:
            try:
                os.chdir(str(label) + '_' + str(option))
            except OSError:
                os.system('mkdir %s' % str(label) + '_' + str(option))
                os.chdir(str(label) + '_' + str(option))
            try:
                os.system("cp ../POSCAR .")
                os.system("cp ../POTCAR .")
            except:
                print 'File system error'
            self.AlterFile("POSCAR",key,option)
            self.DescendPoscar(copy.deepcopy(varmap))
            os.chdir('..')
        try:
            os.system('rm POSCAR')
            os.system('rm POTCAR')
        except:
            print "File system error"
            
    def DescendKpoints(self,varmap):
        varmap = copy.deepcopy(varmap)        
        if not (os.path.isfile('KPOINTS')):
            try:
                os.system('cp %s %s' % (self.Kpoints,'KPOINTS' ))
            except:
                print "File system error"                
        if len(varmap) == 0:
            self.DescendIncar(self.IncarVars)
            return            
        key = varmap.items()[0][0]
        label = key[1:]
        varlist = varmap[key]
        del varmap[key]        
        for option in varlist:
            try:
                os.chdir(str(label) + '_' + str(option))
            except OSError:
                os.system('mkdir %s' % str(label) + '_' + str(option))
                os.chdir(str(label) + '_' + str(option))
            try:
                os.system("cp ../KPOINTS .")
                os.system("cp ../POSCAR .")
                os.system("cp ../POTCAR .")
            except:
                print 'File system error'
            self.AlterFile("KPOINTS",key,option)
            self.DescendKpoints(copy.deepcopy(varmap))
            os.chdir('..')
        try:
            os.system('rm KPOINTS')
            os.system('rm POSCAR')
            os.system('rm POTCAR')
        except:
            print "File system error"
                        
    def DescendIncar(self,varmap):
        varmap = copy.deepcopy(varmap)        
        if not (os.path.isfile('INCAR')):
            try:
                os.system('cp %s %s' % (self.Incar,'INCAR' ))
            except:
                print "File system error"                
        if len(varmap) == 0:
            if self.RunName in os.listdir(os.getcwd()):
                print "Already found " + self.RunName + " in folder " + os.getcwd()
                return
            self.TotalRuns += 1
            os.system('mkdir %s' % self.RunName)
            os.system('mv INCAR %s' % (self.RunName + '/INCAR' ))
            os.system('mv KPOINTS %s' % (self.RunName + '/KPOINTS'))
            os.system('mv POSCAR %s' % (self.RunName + '/POSCAR'))
            os.system('mv POTCAR %s' % (self.RunName + '/POTCAR'))
            return            
        key = varmap.items()[0][0]
        label = key[1:]
        varlist = varmap[key]
        del varmap[key]         
        for option in varlist:
            try: #go to or make dir
                os.chdir(str(label) + '_' + str(option))
            except OSError:
                os.system('mkdir %s' % str(label) + '_' + str(option))
                os.chdir(str(label) + '_' + str(option))
            try: #copy input files from previous level
                os.system("cp ../INCAR .")
                os.system("cp ../KPOINTS .")
                os.system("cp ../POSCAR .")
                os.system("cp ../POTCAR .")
            except:
                print 'File system error'
            self.AlterFile("INCAR",key,option)
            self.DescendIncar(copy.deepcopy(varmap))
            os.chdir('..')
        try:
            os.system('rm INCAR')
            os.system('rm KPOINTS')
            os.system('rm POSCAR')
            os.system('rm POTCAR')
        except:
            print "File system error"

    def AlterFile(self,filepath,frommatch,tomatch):    
        import re
        try:
            fileread = open(filepath,'r')
            text = fileread.read()
            fileread.close()
            filewrite = open(filepath,'w')
            regex = re.compile(str(frommatch))
            text = regex.sub(str(tomatch),text)
            filewrite.write(text)
            filewrite.close()
        except:
            print "There was an error in processing file " + filepath
            print "This is likely due to " + frommatch + " tag not being present in the file."
            print tomatch

    def AlterPotcar(self,element,tag):
        curdir = os.getcwd()
        dir2 = self.PotcarDir+'/'+element+'/'
        os.chdir(dir2)        
        file=open("POTCAR",'r')
        file2=file.readlines()
        file.close()
        file2string=""
        for string in file2:
            file2string=file2string+string
        os.chdir(curdir)
        import re
        try:
            fileread = open('POTCAR','r')
            text = fileread.read()
            fileread.close()
            filewrite = open('POTCAR','w')
            regex = re.compile(str(tag))
            text = regex.sub(str(file2string),text)
            filewrite.write(text)
            filewrite.close()
        except:            
            print "There was an error in processing file " + dir2+"POTCAR"+'/'
            print "This is likely due to " + frommatch + " tag not being present in the file."
            print tomatch

    def FindFolders(self,folder): #Finds all subdirectories
        files = os.listdir(folder)
        for path in files:
            if os.path.isdir(folder+path+"/"):
                self.FolderList.append(folder+path+"/")
                self.FindFolders(folder+path+"/")
    
    def CheckFolders(self):
        for path in self.ToCheckList:
            if path.split("/")[-2] == self.RunName:
                self.CheckedList.append(path)
                print path
                if os.path.exists(path + 'OUTCAR') and finishCheck(path) and self.StepsLessThanNSW(path): 
                    print ('Will skip (finished):'+path)    
                else:
                    self.ToRunList.append(path)             #run only unfinished ones
    def CheckFoldersDOS(self):
        for path in self.ToCheckList:
            if path.split("/")[-2] == self.RunName:
                self.CheckedList.append(path)
                print path
                if os.path.exists(path + 'OUTCAR') and finishCheck(path): 
                    print ('Will skip (finished):'+path)    
                else:
                    self.ToRunList.append(path)  
                
    def CheckForNewRun(self):
        for path in self.CheckedList:
            parpath =  os.path.abspath(os.path.join(path, os.path.pardir))
            if os.path.exists(os.path.join(parpath,newRun)):
                print os.path.join(parpath,newRun) + " already exists."
                if copyFiles:
                    print "Copying " + newRunFile +" from path to current directory."
                    newPath = os.path.join(parpath,newRun) + "/" + newRunFile
                    array = parpath.split("/")
                    newFileName = array[len(array)-2]+array[len(array)-1]+".dat"
                    shutil.copyfile(newPath,mainDir+newFileName)
            else:
                self.ToRunList.append(parpath+"/")

    def finishCheck(self,folder):
        """Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR."""
        lastfolder = os.getcwd()
        os.chdir(folder)
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)    
        return newstring[0].find('Voluntary') > -1 #True/False
    
    def StepsLessThanNSW(self,folder):
        '''Get number of steps in relaxation, and check vs NSW'''
        #get number of steps
        lastfolder = os.getcwd()
        os.chdir(folder)
        if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0: 
            os.chdir(lastfolder) 
            steps = -9999
        oszicar = open('OSZICAR','r')
        laststep = oszicar.readlines()[-1].split()[0]
        oszicar.close()
        os.chdir(lastfolder)  
        try:
            value = int(laststep)
            steps =  value
        except:
            steps =  9999
        #check vs NSW
        proc = subprocess.Popen(['grep','-i','NSW',self.CheckedList[0]+'/INCAR'],stdout=subprocess.PIPE)
        print proc.communicate()[0].split('=')
        NSW = int(proc.communicate()[0].split('=')[-1])
        return steps < NSW # True/False

    def stretchCsPOSCAR(self,runlist):
       '''Takes a POSCAR (in CONTCAR direct coordinates), and moves the C's (first 2 atoms) downward by 
       0.35 of z lattice vector (about 7 ang).  The 2 C's are in lines 9,10, or indices 8,9''' 
       for path in runlist:
           file=open(path + 'POSCAR','r')
           lines=file.readlines()
           file.close()
           if len(lines) == 0:
               return
           c1pos = lines[8].split()
           c2pos = lines[9].split()
           z1 = float(c1pos[2])
           z2 = float(c2pos[2])
           if z1 > 0.5:
               z1 = z1 - 1.0
           if z2 > 0.5:
               z2 = z1 - 1.0
           c1pos = c1pos[0:2]+[str(z1-0.35)]
           c2pos = c2pos[0:2]+[str(z2-0.35)] 
           c1str = c1pos[0] + '  ' + c1pos[1] + '  ' + c1pos[2] + '\n'
           c2str = c2pos[0] + '  ' + c2pos[1] + '  ' + c2pos[2] + '\n'
           lines2 = lines[0:8]+[c1str,c2str]+lines[10:]
           for line in lines2:
               print line
           
           file=open(path + 'POSCAR','w')
#           for i, linei in enumerate(lines2):
           file.writelines(lines2)
           file.close()
       return             
                                              
#    def removeCsPOSCAR(self,runlist):
#       '''Takes a POSCAR (in CONTCAR direct coordinates), and removes the C's (first 2 atoms) so that an 
#       isolated adatomlayer can be run.  The 2 C's are in lines 9,10, or indices 8,9''' 
#       for path in runlist:
#           file=open(path + 'POSCAR','r')
#           lines=file.readlines()
#           file.close()
#           if len(lines) == 0:
#               return
#           c1pos = lines[8].split()
#           c2pos = lines[9].split()
#           z1 = float(c1pos[2])
#           z2 = float(c2pos[2])
#           if z1 > 0.5:
#               z1 = z1 - 1.0
#           if z2 > 0.5:
#               z2 = z1 - 1.0
#           c1pos = c1pos[0:2]+[str(z1-0.35)]
#           c2pos = c2pos[0:2]+[str(z2-0.35)] 
#           c1str = c1pos[0] + '  ' + c1pos[1] + '  ' + c1pos[2] + '\n'
#           c2str = c2pos[0] + '  ' + c2pos[1] + '  ' + c2pos[2] + '\n'
#           lines2 = lines[0:8]+[c1str,c2str]+lines[10:]
#           for line in lines2:
#               print line
#           
#           file=open(path + 'POSCAR','w')
##           for i, linei in enumerate(lines2):
#           file.writelines(lines2)
#           file.close()
#       return             
     
                                    
           
              
               
        
        
    def addHToPOSCAR(self,path,HAdDist):
       '''Takes a POSCAR (in CONTCAR direct coordinates), and adds an H above the adatom (last line)'''
       file=open(path + 'POSCAR','r')
       lines=file.readlines()
       file.close()
       if len(lines) == 0:
       	   return
#       #add H to element list  #only needed if H is not already there
       try:
           elemsString = lines[5]
       except:
       	return
#       len1 = len(elemsString)
#       elem = elemsString.split()[-1]
#       lenelem = len(elem) # last element (adatom) string
#       newstr = elemsString[:len1-lenelem-2]+'H   '+ elem +'\n'#adds H in element list
#       lines[5] = newstr
#       #add to numbers of atoms
       lines[6] ='    2     2     1\n'  # 2 replaces 1 in 2nd (H) spot, in 6th line
       #add position of new H (must scale position to direct coordinates)
       repeatz = float(lines[4].split()[-1])
       newHz = str(float(lines[11].split()[2])+HAdDist/repeatz)
       newHline = '  '+ lines[11].split()[0]+'  '+lines[11].split()[1] +'  ' + newHz+'\n'
       list1 = lines[0:12] #includes up to 11
       list1.insert(11,newHline) #inserts after 11
       #replace old POSCAR
       file=open(path + 'POSCAR','w')
       for i, linei in enumerate(list1):
           file.writelines(linei)
       file.close()
       return                
    
    def replaceParamIncar(self,pathlist, param, value):
        """Replaces a parameter value in Incar with a new value"""
        lastfolder = os.getcwd()
        for path in pathlist:
                os.chdir(path)
                incarfile = open('INCAR', 'r')
                lines = incarfile.readlines()
                incarfile.close()                
                incarfile = open('INCAR', 'w')                
                for line in lines:
                    line2 = line         
                    if param in line2:
                        line2 = param + ' = ' + value + '\n'
                    incarfile.writelines(line2)
                incarfile.close()
        os.chdir(lastfolder)    
        return
