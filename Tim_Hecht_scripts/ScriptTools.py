import os
import copy
class VaspTools:

	def __init__(self,mainFolder,runname,runtype,poscar,kpoints,incar,potcar,
		poscarvar,kpointsvar,incarvar,potcarvar,potcardir):
		
		self._MainFolder = mainFolder
		self._RunName = runname
		self._RunType = runtype
		self._Poscar = poscar
		self._PoscarVars = poscarvar
		self._Kpoints = kpoints
		self._KpointsVars = kpointsvar
		self._Incar = incar
		self._IncarVars = incarvar
		self._Potcar = potcar
		self._PotcarVars = potcarvar
		self._PotcarDir = potcardir
		
		self._FolderList = []

		self._TotalRuns = 0
	
	
	def BuildNewRun(self):
		os.chdir(self._MainFolder)
		
		for folder in self._RunType:
			try:
				os.chdir(folder)
			except OSError:
				os.system('mkdir %s' % folder)
				os.chdir(folder)
		self._DescendPotcar(self._PotcarVars)
		print "Total Number of runs created:"
		print self._TotalRuns
		
				
	def Run(self):
		pass
		
	def _DescendPotcar(self,varmap):
		varmap = copy.deepcopy(varmap)
		
		if not (os.path.isfile('POTCAR')):
			try:
				os.system('cp %s %s' % (self._Potcar,'POTCAR' ))
			except:
				print "File system error"
				
		if len(varmap) == 0:
			self._DescendPoscar(self._PoscarVars)
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
				os.system("cp ../POTCAR .")
			except:
				print 'File system error'
			self._AlterPotcar(option,key)
			self._DescendPotcar(copy.deepcopy(varmap))
			os.chdir('..')
		try:
			os.system('rm POTCAR')
		except:
			print "File system error"
				
	def _DescendPoscar(self,varmap):
		varmap = copy.deepcopy(varmap)
		
		if not (os.path.isfile('POSCAR')):
			try:
				os.system('cp %s %s' % (self._Poscar,'POSCAR' ))
			except:
				print "File system error"
				
		if len(varmap) == 0:
			self._DescendKpoints(self._KpointsVars)
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
			self._AlterFile("POSCAR",key,option)
			self._DescendPoscar(copy.deepcopy(varmap))
			os.chdir('..')
		try:
			os.system('rm POSCAR')
			os.system('rm POTCAR')
		except:
			print "File system error"
			
	def _DescendKpoints(self,varmap):
		varmap = copy.deepcopy(varmap)
		
		if not (os.path.isfile('KPOINTS')):
			try:
				os.system('cp %s %s' % (self._Kpoints,'KPOINTS' ))
			except:
				print "File system error"
				
		if len(varmap) == 0:
			self._DescendIncar(self._IncarVars)
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
			self._AlterFile("KPOINTS",key,option)
			self._DescendKpoints(copy.deepcopy(varmap))
			os.chdir('..')
		try:
			os.system('rm KPOINTS')
			os.system('rm POSCAR')
			os.system('rm POTCAR')
		except:
			print "File system error"

			
	def _DescendIncar(self,varmap):
		varmap = copy.deepcopy(varmap)
		
		if not (os.path.isfile('INCAR')):
			try:
				os.system('cp %s %s' % (self._Incar,'INCAR' ))
			except:
				print "File system error"
				
		if len(varmap) == 0:
			if self._RunName in os.listdir(os.getcwd()):
				print "Already found " + self._RunName + " in folder " + os.getcwd()
				return
			self._TotalRuns += 1
			os.system('mkdir %s' % self._RunName)
			os.system('mv INCAR %s' % (self._RunName + '/INCAR' ))
			os.system('mv KPOINTS %s' % (self._RunName + '/KPOINTS'))
			os.system('mv POSCAR %s' % (self._RunName + '/POSCAR'))
			os.system('mv POTCAR %s' % (self._RunName + '/POTCAR'))
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
				os.system("cp ../INCAR .")
				os.system("cp ../KPOINTS .")
				os.system("cp ../POSCAR .")
				os.system("cp ../POTCAR .")
			except:
				print 'File system error'
			self._AlterFile("INCAR",key,option)
			self._DescendIncar(copy.deepcopy(varmap))
			os.chdir('..')
		try:
			os.system('rm INCAR')
			os.system('rm KPOINTS')
			os.system('rm POSCAR')
			os.system('rm POTCAR')
		except:
			print "File system error"
			
			
	def _AlterFile(self,filepath,frommatch,tomatch):	
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

	def _AlterPotcar(self,element,tag):
		curdir = os.getcwd()
		os.chdir(self._PotcarDir)
    		os.chdir(element)
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
			
		
	def _FindFolders(self,folder):
		files = os.listdir(folder)
		for path in files:
			if os.path.isdir(folder+path+"/"):
				self._FolderList.append(folder+path+"/")
				self._FindFolders(folder+path+"/")
	
		
class ContinuationRunSet:
	def __init__(self):
		print("Hello")
