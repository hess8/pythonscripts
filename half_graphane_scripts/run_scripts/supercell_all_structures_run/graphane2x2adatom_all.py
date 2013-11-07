# Generates all structures possible in a 2x2 graphane supercell, with adatoms replacing H 
# Creates all vasp input files, and perhaps aflow.in
# 
from numpy import array, arccos, dot, pi, zeros, floor, sum
from numpy.linalg import norm
from random import random
import os, subprocess, sys, time
from runtools import writePoscar, writejobfile

#maindir = '/fslhome/bch/cluster_expansion/hexagonal/aflow2x2adatoms/'
maindir = '/fslhome/bch/cluster_expansion/hexagonal/2x2adatoms/'
aflow = False #1 to turn aflow prep on or off
#run = 'relaxfinal'
run = 'relax'

# L = input('\nBond length? ')
Nsuper1 = 2 # N of supercells in two directions
Nsuper2 = 2
NCcell = 2 # number of C atoms in primitive unit cell
NHcell = 2 # starting H locations in unit cell
atomslist = ['C','H','W']
dAd = 2.2  # Adatom  distance from plane

nsites = NCcell*(Nsuper1+Nsuper2)
typetop = zeros((8))
nstructs = 2**nsites
for istruct in range(nstructs): # 
    typetopstr = bin(istruct)[2:] # binary: 0 if H, 1 if adatom 
    print istruct
    for i in range(len(typetopstr)) :
       typetop[i+nsites-len(typetopstr)] = int(typetopstr[i])#fill out all spaces in the binary representation
    print  typetop #binary representation of H/adatoms   
    dir = maindir + 'struct%s/' %istruct 
    if not os.path.isdir(dir):
        os.mkdir(dir)       
    if not os.path.isdir(dir+run):
        os.mkdir(dir+run)     
    os.chdir(dir+run)
    if not os.path.exists('converged.dat'): # do only dirs that have not converged!!!!       
        os.system('rm slurm*.out')
        if run == 'relax' and (not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0): 
            # first run, write POSCAR
            listused = getAtomsUsed(atomslist,istruct,nstructs)         
            writePoscar(typetop,listused,istruct, nstructs,Nsuper1,Nsuper2,NCcell)
            writePotcar(listused)
            #create POTCAR          
        elif run == 'relax' and os.path.exists('CONTCAR'):
            os.system('cp POSCAR POSCAR%s' % time.strftime("%X")) #appends 24hr time
            os.system('cp CONTCAR POSCAR')     
        elif run == 'relaxfinal':
            os.system('cp ../relax/CONTCAR POSCAR')
            os.system('cp ../relax/POSCAR POSCAR.orig')        
            os.system('cp ../relax/POTCAR .')                         
        if aflow:
            os.system('aconvasp --poscar2aflowin < POSCAR > aflow.in')
        else:  #Copy vasp input files              
            os.system('cp ../../../vaspinput/%s/KPOINTS .' % run)
            os.system('cp ../../../vaspinput/%s/INCAR .' % run)
            os.system('cp ../../../vaspinput/%s/job .' % run)
        writejobfile('struct'+str(istruct),atomslist[-1]+run,'vasp533')
        os.system('sbatch job')
        print 'Submitted job for struct%s' % str(istruct)
#create aflow.in

os.chdir(maindir)
#os.system("find `pwd` -name 'aflow.in' > jobs2run")

print 'done'

