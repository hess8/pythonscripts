#analysis script
mainDir = '/bluehome/bch/vasprun/graphene.structures/half_graphane/'
subdir = 'final_relax'
dir = mainDir + subdir + '/'
isolateddir = mainDir + 'isolated/'
initialdir = mainDir + 'initial_relax/'
import sys


file = open('outputatoms.dat','r')
atomsfile = file.readlines()
atomsdata = {}
file.close()

file = open('output1elementhollow.dat','r')
sysfile = file.readlines()
sysdata={}
file.close()

file = open('outputstretch.dat','r')
stretchfile = file.readlines()
stretchdata={}
file.close()

file = open('outputconvergence.dat','r')
convergencefile = file.readlines()
convergencedata={}
file.close

file = open('outputdistance.dat','r')
distancefile = file.readlines()
distancedata={}
file.close

file = open('outputcarbon.dat','r')
carbonfile = file.readlines()
carbondata={}
file.close

file = open('outputexpressure.dat','r')
expressurefile = file.readlines()
expressuredata={}
file.close

file = open('outputexfoliation.dat','r')
exfoliationfile = file.readlines()
exfoliationdata={}
file.close

file = open('outputlatticevector.dat','r')
latticevectorfile = file.readlines()
latticevectordata={}
file.close

for i in range(len(atomsfile)):
  atomsdata[atomsfile[i].split()[0]]=float(atomsfile[i].split()[1])

for i in range(len(sysfile)):
  try:
      sysdata[sysfile[i].split()[0]]=float(sysfile[i].split()[1])
  except:
      sysdata[sysfile[i].split()[0]]=None
  
for i in range(len(stretchfile)):
  try:
    stretchdata[stretchfile[i].split()[0]]=float(stretchfile[i].split()[1])
  except:
    stretchdata[stretchfile[i].split()[0]]=None


for i in range(len(convergencefile)):
  try:
    convergencedata[convergencefile[i].split()[0]]=float(convergencefile[i].split()[1])
  except:
    convergencedata[convergencefile[i].split()[0]]=None

for i in range(len(distancefile)):
  try:
    distancedata[distancefile[i].split()[0]]=float(distancefile[i].split()[1])
  except:
    distancedata[distancefile[i].split()[0]]=None

for i in range(len(carbonfile)):
  try:
    carbondata[carbonfile[i].split()[0]]=float(carbonfile[i].split()[1])
  except:
    carbondata[carbonfile[i].split()[0]]=None

for i in range(len(expressurefile)):
  try:
    expressuredata[expressurefile[i].split()[0]]=float(expressurefile[i].split()[1])
  except:
    expressuredata[expressurefile[i].split()[0]]=None

for i in range(len(exfoliationfile)):
  try:
    exfoliationdata[exfoliationfile[i].split()[0]]=float(exfoliationfile[i].split()[1])
  except:
    exfoliationdata[exfoliationfile[i].split()[0]]=None
    
for i in range(len(latticevectorfile)):
  try:
    latticevectordata[latticevectorfile[i].split()[0]]=float(latticevectorfile[i].split()[1])
  except:
    latticevectordata[latticevectorfile[i].split()[0]]=None

file = open('outputanalysishollow.csv','w')
file.write('System,Number of Atoms,Calculated Energy,Calulated Energy (Stretched),Adatom Monolayer Energy,Lattice Vector of Monolayer,Zero Energy,Zero Energy with Graphene,Atomization Energy,Atomization Energy per Atom,Binding Energy,BE per Adatom,Interplane Binding Energy,IBE Per Adatom,Exfoliation Energy,EE Per Adatom,Largest Force on Atom,Distance Between Carbon and Adatom,Carbon Distance out of Plane,External Pressure\n')
for label in sysdata.keys():
  templist =  label.split(":")
  numAtoms = 2
  totEnergy = 2*atomsdata["C"]
  bindEnergy = -18.44
  for item in templist[0].split("-"):
    totEnergy = totEnergy + atomsdata[item]
    bindEnergy=bindEnergy+atomsdata[item]
    numAtoms = numAtoms+1
  for item in templist[1].split("-"):
    if item == "":
      continue
    totEnergy = totEnergy + atomsdata[item]
    bindEnergy = bindEnergy + atomsdata[item]
    numAtoms = numAtoms+1
    
  try:
    string=label + "," #system
    string += str(numAtoms)+"," #numatoms
    string += str(sysdata[label])+"," #calculated Energy
    string += str(stretchdata[label])+"," #stretch energy
    if numAtoms == 3:
        string += str(exfoliationdata[label])+"," #monolayer energy
        string += str(latticevectordata[label])+"," #Lattice Vector
    else:
        string += str(2*exfoliationdata[label.split(":")[0]+":"])+"," #monolayer energy
        string += str(latticevectordata[label.split(":")[0]+":"])+"," #Lattice Vector
    string += str(totEnergy)+"," #zeroenergy
    string += str(bindEnergy)+","#zero energy with graphene as zero
    string += str(sysdata[label]-totEnergy)+"," #Atomization energy
    string += str((sysdata[label]-totEnergy)/float(numAtoms))+"," #AE per atom
    string += str(sysdata[label]-bindEnergy)+"," #binding energy
    string += str((sysdata[label]-bindEnergy)/float(numAtoms-2))+"," #BE per adatom
    string += str(sysdata[label]-stretchdata[label])+"," #IBE
    string += str((sysdata[label]-stretchdata[label])/(numAtoms-2))+"," #IBE per adatom
    if numAtoms == 3:
        string += str(-18.44+exfoliationdata[label.split(":")[0]+":"]-sysdata[label])+","
        string += str(-18.44+exfoliationdata[label.split(":")[0]+":"]-sysdata[label])+","
    else:
        string += str(-18.44+2*exfoliationdata[label.split(":")[0]+":"]-sysdata[label])+","
        string += str((-18.44+2*exfoliationdata[label.split(":")[0]+":"]-sysdata[label])/2)+","
    string += str(convergencedata[label])+"," #convergence
    string += str(distancedata[label])+","
    string += str(carbondata[label])+","
    string += str(expressuredata[label])+"\n"
  except:
    string = label + " Error\n"
  file.write(string)


file.close()
