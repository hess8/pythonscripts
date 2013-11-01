# Generates all structures possible in a 2x2 graphane supercell, with adatoms replacing H 
# 
from numpy import array, arccos, dot, pi, zeros, floor, sum
from numpy.linalg import norm
from random import random
import os, subprocess, sys

dir = '/fslhome/bch/cluster_expansion/hexagonal/2x2adatoms'
os.chdir(dir)

# L = input('\nBond length? ')
Nsuper1 = 2 # N of supercells in two directions
Nsuper2 = 2
NCcell = 2 # number of C atoms in unit cell
NCall = 2*Nsuper1*Nsuper2
NHcell = 2 # starting H locations in unit cell
atomslist = 'C H W \n'
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
rc1 =      array([1.46545,      0.00000,      0.22856])        
rc2 =      array([2.93090,      0.00000,     -0.22856])         
rh1 =      array([1.46545,      0.00000,      1.33850])      
rh2 =      array([2.93090,      0.00000,     -1.33850 ])      

a1s = Nsuper1*a1new # superlattice vectors, including break
a2s =  Nsuper2*a2new  #  
r = zeros((Nsuper1,Nsuper2,NCcell,3))
rh = zeros((Nsuper1,Nsuper2,NCcell,3))

# C 
for j in range(Nsuper1):
    for k in range(Nsuper2):
        r[j,k,0,:] = rc1 + (j-1)*a1new + (k-1)*a2new
        r[j,k,1,:] = rc2 + (j-1)*a1new + (k-1)*a2new

#  Hydrogen sites for Chair graphene
for j in range(Nsuper1):
    for k in range(Nsuper2):      
        rh[j,k,0,:] = rh1 + (j-1)*a1new + (k-1)*a2new
        rh[j,k,1,:] = rh2 + (j-1)*a1new + (k-1)*a2new



typetop = zeros((8))
for istruct in range(2**(2*(Nsuper1+Nsuper2))): # in each loop we add one more adatom to the structure, sequentially
    typetopstr = bin(istruct)[2:] # binary: 0 if H, 1 if adatom 
    print typetopstr
    for i in range(len(typetopstr)) :
       typetop[i] = int(typetopstr[i])
    print  typetop
    NAd = sum(typetop)
    NHall = NCall - NAd #rest of top sites covered by H       
    rad = zeros((NAd,3))
    #Create cell positions
       
    pscr=open('POSCAR'+str(istruct) ,'w')
    pscr.write(atomslist)
    pscr.write('1.0 \n')
    #  lattice vectors
    
    pscr.write('%12.8f %12.8f %12.8f  \n' % (a1s[0], a1s[1], a1s[2]))
    pscr.write('%12.8f %12.8f %12.8f  \n' % (a2s[0], a2s[1], a2s[2]))
    pscr.write('%12.8f %12.8f %12.8f  \n' % (a3s[0], a3s[1], a3s[2]))
    
    pscr.write('%g %g %g\n' % (NCall, NHall, NAd))
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
                   print whichsite
                   if typetop[whichsite] == 1: #this is adatom
                       print 'Site %s is adatom' % str(whichsite)
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
print 'done'

