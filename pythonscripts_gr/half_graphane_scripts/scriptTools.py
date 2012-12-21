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
            print 'path', path
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
        print key
        label = key[1:] #keyword @adatom        
        varlist = varmap[key]
        del varmap[key]
        print key
        for option in varlist:
            print (option)
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
        os.chdir(self.PotcarDir+'/'+element)        
        # don't know why these are indented; error without it
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
            print "There was an error in processing file " + filepath
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
#            print("CHECK NEXT LINE")
#            print(path.split("/")[-2])
            if path.split("/")[-2] == self.RunName:
                self.CheckedList.append(path)
                if os.path.exists(path + 'OUTCAR') and self.FinishCheck(path): 
                    print ('Will skip (finished):'+path)    
                else:
                    self.ToRunList.append(path)             #run only unfinished ones

                
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

    def FinishCheck(self,folder):
    #        """Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR."""
        lastfolder = os.getcwd()
        os.chdir(folder)
        print folder
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        print newstring[0]
        os.chdir(lastfolder)    
        return newstring[0].find('Voluntary') > -1 #True/False

    def Contcar2Poscar(self,bestfolder,folder):
        lastfolder = os.getcwd()
        os.system('cp ', bestfolder + 'CONTCAR POSCAR')
        lastfolder = os.getcwd()
        os.chdir(folder)
        os.system('cp ', bestfolder + 'CONTCAR POSCAR')
        os.chdir(lastfolder)
