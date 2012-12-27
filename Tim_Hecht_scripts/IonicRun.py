#Specify Directory to use
mainDir = "/bluehome/thecht/TransitionMetals/"

#Specify Potcar Directory
potcardir = "/bluehome/thecht/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/"

#Specify the type of run
runType = ['4x4','Hollow','Single']

#Specify the name of the type of run
runName = "IonicRun"

#Specify a Poscar file
poscar = mainDir + 'Poscar/4x4HollowSiteSingle.poscar'

#Specify a KPoints file
kpoints = mainDir + 'Kpoints/grapheneionicrun.kpoints'

#Specify an Incar file
incar = mainDir + 'Incar/grapheneionicrun.incar'

#Specify a Potcar file
potcar = mainDir + 'Potcar/CH.potcar'

#Specify Poscar variables
poscarVariables = {
'@distance':
	[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0]
}

#Specify KPoints variables
kpointVariables = {

}

#Specify Incar variables
incarVariables = {
'@ISPIN':
	['1','2'],
'@IBRION':
	['-1','2'],
'@ISIF':
	['3','4']

}

#Specify Potcar Elements
elementList = {
'@Adatom':
	['Hf','La']
}

import ScriptTools

tools = ScriptTools.VaspTools(mainDir,runName,runType,poscar,kpoints,incar,
	potcar,poscarVariables,kpointVariables,incarVariables,elementList,potcardir)

tools.BuildNewRun()

print "Done"
