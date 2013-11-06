

def writejobfile(struct,nameadd,vaspexec):
    '''read from a template one level up from maindir, and put dir in job name ('myjob' is replaced).  
    Choose executable label (should be linked in bin file)'''
    print os.getcwd()
    file1 = open('job','r')
    jobfile = file1.readlines()
    file1.close
    for i in range(len(jobfile)):
        jobfile[i]=jobfile[i].replace('myjob', struct+nameadd)
        if '>' in jobfile[i]:
            jobfile[i]= vaspexec + ' > out'
    file2 = open('job','w')
    file2.writelines(jobfile) 
    file2.close()
    return 
# Generates all structures possible in a 2x2 graphane supercell, with adatoms replacing H 
# Creates all vasp input files, and aflow.in
# 
from numpy import array, arccos, dot, pi, zeros, floor, sum
from numpy.linalg import norm
from random import random
import os, subprocess, sys

#maindir = '/fslhome/bch/cluster_expansion/hexagonal/aflow2x2adatoms/'
maindir = '/fslhome/bch/cluster_expansion/hexagonal/2x2adatoms/'
aflow = False #1 to turn aflow prep on or off
#run = 'relaxfinal'
run = 'relax'

# L = input('\nBond length? ')
Nsuper1 = 2 # N of supercells in two directions
Nsuper2 = 2
NCcell = 2 # number of C atoms in unit cell
NCall = 2*Nsuper1*Nsuper2
NHcell = 2 # starting H locations in unit cell
atomslist = ['C','H','W']
dAd = 2.2  # Adatom  distance from plane
                                                                                                                          
#Lattice Vectors from vasp primitive cell relaxation 
a1 =    array([2.13129, -1.2305, 0 ])     
a2 =    array([2.13129, 1.2305, 0  ])    
a3s =    array([.00000000,  0.000000000, 15.00]) 
#this one will be unchanged in supercell
theta = arccos(dot(a1,a2)/norm(a1)/norm(a2)) #in case of distorted hexagonal lattice
thetadeg = theta*180/pi
# Doubling size of unit cell to allow boat structures. New primitive cell can be made of a diagonal and a horizontal vector
# a1new = a1 + a2#5# horizontal...need to double number of atoms in unit cell.
# a2new = a2 - a1 #(gives a vertical lattice vector)
a1new = a1 # Don't change
a2new = a2 # 
#  Atom positions from vasp primitive cell relaxation
#for "graphane" start of plane (diamond like)
#rc1 =      array([1.46545,      0.00000,      0.22856])        
#rc2 =      array([2.93090,      0.00000,     -0.22856])         
#rh1 =      array([1.46545,      0.00000,      1.33850])      
#rh2 =      array([2.93090,      0.00000,     -1.33850 ])    

#For "flat start of plane"
rc1 =      array([1.46545,      0.00000,      0.0])        
rc2 =      array([2.93090,      0.00000,     -0.0])         
rh1 =      array([1.46545,      0.00000,      1.1])      
rh2 =      array([2.93090,      0.00000,     -1.1 ])   

a1s = Nsuper1*a1new # superlattice vectors, including break
a2s =  Nsuper2*a2new  #  
r = zeros((Nsuper1,Nsuper2,NCcell,3))
rh = zeros((Nsuper1,Nsuper2,NCcell,3))

# C 
for j in range(Nsuper1):
    for k in range(Nsuper2):
        r[j,k,0,:] = rc1 + (j)*a1new + (k)*a2new
        r[j,k,1,:] = rc2 + (j)*a1new + (k)*a2new

#  Hydrogen sites for Chair graphene
for j in range(Nsuper1):
    for k in range(Nsuper2):      
        rh[j,k,0,:] = rh1 + (j)*a1new + (k)*a2new
        rh[j,k,1,:] = rh2 + (j)*a1new + (k)*a2new

nsites = NCcell*(Nsuper1+Nsuper2)
typetop = zeros((8))
nstructs = 2**nsites
for istruct in range(nstructs): # in each loop we add one more adatom to the structure, sequentially
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
        if run == 'relax' and (not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0): # first run, write POSCAR
            NAd = sum(typetop)
            NHall = NCall - NAd #rest of top sites covered by H       
            rad = zeros((NAd,3))
            #Create cell positions
            pscr=open('POSCAR','w')
            if istruct ==0:
                listused = [atomslist[0], atomslist[1]] # all H
            elif istruct == nstructs-1: #all adatoms
                listused = [atomslist[0], atomslist[2]]
            else:
                listused =atomslist 
            for atom in listused:
                pscr.write('%s ' %atom)
        #    pscr.write('\n')     
            pscr.write('\n1.0\n')
            #  lattice vectors
            
            pscr.write('%12.8f %12.8f %12.8f  \n' % (a1s[0], a1s[1], a1s[2]))
            pscr.write('%12.8f %12.8f %12.8f  \n' % (a2s[0], a2s[1], a2s[2]))
            pscr.write('%12.8f %12.8f %12.8f  \n' % (a3s[0], a3s[1], a3s[2]))
            
            pscr.write('%g ' % NCall)
            if NHall>0:
                pscr.write('%g ' % NHall)
            if NAd>0:
                pscr.write('%g ' % NAd)
            pscr.write('\n')      
            pscr.write('Cartesian \n')
        
            #print (#no atom type numbers for POSCAR...just a list of positions
            for j in range(Nsuper1):
                for k in range(Nsuper2):
                    for iat in range(NCcell):
                        pscr.write('%12.8f %12.8f %12.8f \n' % ( r[j,k,iat,0], r[j,k,iat,1], r[j,k,iat,2]))
                #print hydrogen atoms before adatoms
            adcount = 0
            whichsite = 0
            for j in range(Nsuper1):
                for k in range(Nsuper2):
                    for iat in range(NHcell):                  
                           if typetop[whichsite] == 1: #this is adatom
        #                       print 'Site %s is adatom' % str(whichsite)
                               if rh[j,k,iat,2] > 0: #position is above the plane
                                   rad[adcount,:] = [rh[j,k,iat,0], rh[j,k,iat,1], dAd]
                               else:
                                   rad[adcount,:] = [rh[j,k,iat,0], rh[j,k,iat,1], -dAd]
                               adcount += 1
                           else:
                               pscr.write('%12.8f %12.8f %12.8f \n' % ( rh[j,k,iat,0], rh[j,k,iat,1], rh[j,k,iat,2]))
                           whichsite += 1
            #print adatom positions
            for j in range(adcount):
                pscr.write('%12.8f %12.8f %12.8f \n' % ( rad[j,0], rad[j,1], rad[j,2]))
            pscr.close()
            #create POTCAR       
            potstr = 'cat ' 
            for atom in listused:
                potstr += '~/vaspfiles/src/potpaw_PBE/%s/POTCAR ' %atom
            potstr += '> POTCAR'
            os.system(potstr)   
        elif run == 'relax' and os.path.exists('CONTCAR'):
#            os.system('cp POSCAR POSCAR.orig')
            os.system('cp CONTCAR POSCAR')     
        elif run == 'relaxfinal':
            os.system('cp ../relax/CONTCAR POSCAR')
            os.system('cp ../relax/POSCAR POSCAR.orig')        
            os.system('cp ../relax/POTCAR .')                         
        if aflow:
            dummy=0
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

