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

    def __init__(self,mainDir,inputDir,runname,runtype,inputlist,toCheckList,checkedList, toRunList ):
        self.MainDir = mainDir
        self.InputDir = inputDir
        self.RunName = runname
        self.RunType = runtype
        self.InputList = inputlist       
        self.FolderList = []
        self.DirsCreated = 0
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
        os.chdir(self.MainDir)
        path1 = self.MainDir+self.RunType+'/'
        try:
            os.chdir(path1)
        except:
            os.system('mkdir %s' % path1)
            os.chdir(path1)
        self.DirLoop(self.InputList) 
        print "Number of run folders created:"
        print self.DirsCreated
        
    def DirLoop(self,inputlist):
        varSets = []
        for list1 in inputlist:
            varSets = varSets + list1  #flatten the inputlist. Each  is of the form  ['@NFITSTRUC','CS.in',[64,128]],  Includes tag and file, vars
        for list1 in varSets:
            print list1
            inputpath = self.InputDir+list1[1] #file path
            os.system('cp '+ inputpath + ' .')         # copy inputfiles to first level
        self.DescendOption(varSets)
        print 'Descend done'
                        
    def DescendOption(self,listRemain): #"outer loop" of file structure. Creates run folder if needed
        listRemain = copy.deepcopy(listRemain)
        print 'Rem:', len(listRemain),listRemain
        if len(listRemain) == 0: #done
            if self.RunName in os.listdir(os.getcwd()):
                print "Already found " + self.RunName + " in folder " + os.getcwd()
                return
            self.DirsCreated += 1
            os.system('mkdir %s' % self.RunName)
            os.system('ls')
            os.system('mv *.in %s/' % self.RunName)                                  
            return               
        tempList = listRemain[0] #for this level of dirs
        tag = tempList[0]
        tagtext= tag.lower()[1:] #drops the @ symbol and uses only lowercase for dir name
        filename = tempList[1]
        filepath = self.InputDir+'/'+filename
        levelVars = tempList[2]
        for var in levelVars:   
            print tagtext, var        
            dirVar = str(tagtext) + '.' + str(var)
            try:
                os.chdir(dirVar)
            except OSError:
                os.system('mkdir %s' % dirVar)
                os.chdir(dirVar)
            try:
                os.system('cp ../*.in .')
            except:
                print 'File system error'
            self.AlterFile(filename,tag,var)
            self.DescendOption(copy.deepcopy(listRemain[1:]) ) #removes the option that is done, then does recursion.
            os.chdir('..')
        try:
                os.system('rm *.in')
        except:
            print "File system error\n"                           

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

#    def AlterCsIn(self,option,tag):
#        curdir = os.getcwd()
#        dir2 = self.CsInFile+'/'+option+'/'
#        os.chdir(dir2)        
#        file=open("CS.in",'r')
#        file2=file.readlines()
#        file.close()
#        file2string=""
#        for string in file2:
#            file2string=file2string+string
#        os.chdir(curdir)
#        import re
#        try:
#            fileread = open('CS.in','r')
#            text = fileread.read()
#            fileread.close()
#            filewrite = open('CS.in','w')
#            regex = re.compile(str(tag))
#            text = regex.sub(str(file2string),text)
#            filewrite.write(text)
#            filewrite.close()
#        except:            
#            print "There was an error in processing file " + dir2+"CS.in"+'/'
#            print "This is likely due to " + frommatch + " tag not being present in the file."
#            print tomatch
#
#    def FindFolders(self,folder): #Finds all subdirectories
#        files = os.listdir(folder)
#        for path in files:
#            if os.path.isdir(folder+path+"/"):
#                self.FolderList.append(folder+path+"/")
#                self.FindFolders(folder+path+"/")
#    
#    def CheckFolders(self):
#        for path in self.ToCheckList:
#            if path.split("/")[-2] == self.RunName:
#                self.CheckedList.append(path)
#                if os.path.exists(path + 'OUTCAR') and self.FinishCheck(path) and self.StepsLessThanNSW(path): 
#                    print ('Will skip (finished):'+path)    
#                else:
#                    self.ToRunList.append(path)             #run only unfinished ones
#                
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
#
#    def replaceParamIncar(self,pathlist, param, value):
#        """Replaces a parameter value in Incar with a new value"""
#        lastfolder = os.getcwd()
#        for path in pathlist:
#                os.chdir(path)
#                incarfile = open('INCAR', 'r')
#                lines = incarfile.readlines()
#                incarfile.close()                
#                incarfile = open('INCAR', 'w')                
#                for line in lines:
#                    line2 = line         
#                    if param in line2:
#                        line2 = param + ' = ' + value + '\n'
#                    incarfile.writelines(line2)
#                incarfile.close()
#        os.chdir(lastfolder)    
#        return
