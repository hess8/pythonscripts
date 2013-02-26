import os, subprocess, copy, string  

def nstrip(list):
    '''Strips off /n'''
    import string
    list2 = []
    for string1 in list:   
        string2 = string1.strip("\n")
        list2.append(string2)
    return list2

class ceTools:

    def __init__(self,mainDir,inputDir,vaspDataDir, runname,runtype,inputlist,toCheckList,checkedList, toRunList ):
        self.MainDir = mainDir
        self.InputDir = inputDir
        self.VaspDataDir = vaspDataDir
        self.RunName = runname
        self.RunType = runtype
        self.InputList = inputlist       
        self.FolderList = []
        self.DirsCreated = 0
        self.ToCheckList = toCheckList
        self.CheckedList = checkedList
        self.ToRunList = toRunList
        self.nstruc = 0
        self.nfit = 0 
        self.n2 = 0 
        self.n3 = 0 
        self.n4 = 0 
        self.n5 = 0 
        self.n6 = 0       
       
    def AddToList(self,folder):
        files = os.listdir(folder)
        for path in files:
            if os.path.isdir(folder+path+"/"):
                self.ToCheckList.append(folder+path+"/")
                self.AddToList(folder+path+"/")        
    
    def BuildNewRun(self):
        os.chdir(self.MainDir)
        path1 = self.MainDir+self.RunType+'/'
        try:
            os.chdir(path1)
        except:
            os.system('mkdir %s' % path1)
            os.chdir(path1)
        varSets = []
        for list1 in self.InputList:
            varSets = varSets + list1  #flatten the inputlist. Each  is of the form  ['@NFITSTRUC','CS.in',[64,128]],  Includes tag and file, vars
        for list1 in varSets:
            inputpath = self.InputDir+list1[1] #file path
            os.system('cp '+ inputpath + ' .')         # copy inputfiles to first level
        self.DescendOption(varSets)
        print "Number of run folders created:", self.DirsCreated
 
    def DescendOption(self,listRemain):
        '''Creates dir structure if needed and copies input files into run folders'''
        if len(listRemain) == 0: #done
            if self.RunName in os.listdir(os.getcwd()):
                print "Already found " + self.RunName + " in folder " + os.getcwd()
                return
            self.DirsCreated += 1
            os.system('mkdir %s' % self.RunName)
            os.system('mv *.in %s/' % self.RunName)
#            print 'copy',  self.VaspDataDir, os.getcwd()
            os.system('cp %s* %s/' % (self.VaspDataDir,self.RunName))    #move structure data from vasp to run folder                                                     
            return               
        tempList = listRemain[0] #for this level of dirs
        tag = tempList[0]
        tagtext= tag.lower()[1:] #drops the @ symbol and uses only lowercase for dir name
        filename = tempList[1]
        filepath = self.InputDir+'/'+filename
        levelVars = tempList[2]
        for var in levelVars:
            if tagtext == 'n2body': #need to seed n2 for other levels if using "grow" method
                self.n2 = var       
            dirVar = str(tagtext) + '_' + str(var)
            try:
                os.chdir(dirVar)
            except OSError:
                os.system('mkdir %s' % dirVar)
                os.chdir(dirVar)
            try:
                os.system('cp ../*.in .')
            except:
                print 'File system error'
            if tag[0] == '@':
                self.AlterFile(filename,tag,var) #alters input files for uncle
            else:
                self.CalculateAndAlter(tag,var)
            self.DescendOption(listRemain[1:]) #removes the option that is done, then does recursion.
            os.chdir('..')
        try:
                os.system('rm *.in')
        except:
            print "File system error\n"                           

    def CalculateAndAlter(self,tag, var):
        '''Calculates values and alters files'''
        if tag == '*grow': # means need to calculate
            for order in [3,4,5,6]:
                temptag = '@N%sBODY' % str(order)
                norder = str(int(self.n2*var**(order-2)))
                self.AlterFile('lat.in',temptag,norder) #alters input files for uncle
                      
    def RemoveBadRuns(self,templist1):
        '''Removes folders if number of structures is greater than number of clusters'''
        templist = templist1[:]
        ndrop = 0
        for path in self.ToRunList:
            self.getRunValues(path)
            if self.nstruc > self.n2+self.n3+self.n4+self.n5+self.n6:
                print 'Removing bad run', path
                templist.remove(path)
                os.system('rm -r %s' %path) 
                ndrop += 1             
#        templist1 = templist[:]
        print '\nRemoved %s folders from list and file structure' % ndrop
        return templist[:]
                           
    def getRunValues(self, path):
        '''gets each tag's value from the path'''
        for segment in path.split('/'):
            testsplit = segment.split('_')
            if 'nfitstruc' in testsplit:
                self.nstruc = int(testsplit[1])
            if 'nfits' in testsplit:
                self.nfits = int(testsplit[1])                
            if 'n2body' in testsplit:
                self.n2 = int(testsplit[1])                        
            if 'n3body' in testsplit:
                self.n3 = int(testsplit[1]) 
            if 'n4body' in testsplit:
                self.n4 = int(testsplit[1]) 
            if 'n5body' in testsplit:
                self.n5 = int(testsplit[1]) 
            if 'n6body' in testsplit:
                self.n6 = int(testsplit[1])
            if 'grow' in testsplit:  
                growVar = testsplit[1]
                self.n3 = int(self.n2 * float(growVar))
                self.n4 = int(self.n3 * float(growVar))                                                                                          
                self.n5 = int(self.n4 * float(growVar))
                self.n6 = int(self.n5 * float(growVar))               
                                                               
    def AlterFile(self,filepath,frommatch,tomatch):    
        '''Substitutes input tags with values'''
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
                if os.path.exists(path + 'complete.txt'): 
                    print ('Will skip (complete):'+path)    
                else:
                    self.ToRunList.append(path)             #run only unfinished ones
                
#    def CheckForNewRun(self):
#        for path in self.CheckedList:
#            parpath =  os.path.abspath(os.path.join(path, os.path.pardir))
#            if os.path.exists(os.path.join(parpath,newRun)):
#                print os.path.join(parpath,newRun) + " already exists."
#                if copyFiles:
#                    print "Copying " + newRunFile +" from path to current directory."
#                    newPath = os.path.join(parpath,newRun) + "/" + newRunFile
#                    array = parpath.split("/")
#                    newFileName = array[len(array)-2]+array[len(array)-1]+".dat"
#                    shutil.copyfile(newPath,mainDir+newFileName)
#            else:
#                self.ToRunList.append(parpath+"/")
#
#    def FinishCheck(self,folder):
#        """Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR."""
#        lastfolder = os.getcwd()
#        os.chdir(folder)
#        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
#        newstring = proc.communicate()
#        os.chdir(lastfolder)    
#        return newstring[0].find('Voluntary') > -1 #True/False
#    
#    def StepsLessThanNSW(self,folder):
#        '''Get number of steps in relaxation, and check vs NSW'''
#        #get number of steps
#        lastfolder = os.getcwd()
#        os.chdir(folder)
#        if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0: 
#            os.chdir(lastfolder) 
#            steps = -9999
#        oszicar = open('OSZICAR','r')
#        laststep = oszicar.readlines()[-1].split()[0]
#        oszicar.close()
#        os.chdir(lastfolder)  
#        try:
#            value = int(laststep)
#            steps =  value
#        except:
#            steps =  9999
#        #check vs NSW
#        proc = subprocess.Popen(['grep','-i','NSW',self.CheckedList[0]+'/INCAR'],stdout=subprocess.PIPE)
#        NSW = int(proc.communicate()[0].split('=')[-1])
#        return steps < NSW # True/False      