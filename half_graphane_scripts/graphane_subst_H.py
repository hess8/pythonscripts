# This version allows addition of hydrogen 
# and an adatom for testing one adatom in a supercell.
from numpy import array, arccos, dot, pi, zeros
from numpy.linalg import norm
# L = input('\nBond length? ')
Nsuper1 = 4 # N of supercells in two directions
Nsuper2 = 2
Ncarbon = 4 # number of C atoms in unit cell
NH = 4
NAd = 1 # N adatoms in the entire supercell

#Lattice Vectors from vasp primitive cell relaxatio   
a1 =    array([2.13129, -1.2305, 0 ])     
a2 =    array([2.13129, 1.2305, 0  ])    
a3s =    array([.00000000,  0.000000000, 20.00]) 
#this one will be unchanged in supercell
theta = arccos(dot(a1,a2)/norm(a1)/norm(a2)) #in case of distorted hexagonal lattice
thetadeg = theta*180/pi
# Doubling size of unit cell to allow boat structures. New primitive cell can be made of a diagonal and a horizontal vector
# a1new = a1 + a2#5# horizontal...need to double number of atoms in unit cell.
# a2new = a2 - a1 #(gives a vertical lattice vector)
a1new = a1 # Don't change
a2new = 2*a2 # Because we have 4 in cell

#  Atom positions from vasp primitive cell relaxation
rc1 =      array([1.46545,      0.00000,      0.22856])        
rc2 =      array([2.93090,      0.00000,     -0.22856])         
rh1 =      array([1.46545,      0.00000,      1.33850])      
rh2 =      array([2.93090,      0.00000,     -1.33851 ])      


dAd = 1.8  # Adatom  distance from plane
                                                                                                                          
a1s = Nsuper1*a1new # superlattice vectors, including break
a2s =  Nsuper2*a2new  #  
r = zeros((Nsuper1,Nsuper2,4,3))
rh = zeros((Nsuper1,Nsuper2,4,3))

#Create cell positions
   
file1=open('POSCAR' ,'w')
file1.write('C H Ti \n') 
file1.write('1 \n')
#  lattice vectors

file1.write('%8.6f %8.6f %8.6f  \n' % (a1s[0], a1s[1], a1s[2]))
file1.write('%8.6f %8.6f %8.6f  \n' % (a2s[0], a2s[1], a2s[2]))
file1.write('%8.6f %8.6f %8.6f  \n' % (a3s[0], a3s[1], a3s[2]))

file1.write('%g %g %g\n' % (Ncarbon*Nsuper1*Nsuper2, NH*Nsuper1*Nsuper2, NAd))
file1.write('Cartesian \n')
# C 
for j in range(Nsuper1):
    for k in range(Nsuper2):
        r[j,k,0,:] = rc1 + (j-1)*a1new + (k-1)*a2new
        r[j,k,1,:] = rc2 + (j-1)*a1new + (k-1)*a2new
        r[j,k,2,:] = (r[j,k,0,:]) + a2
        r[j,k,3,:] = (r[j,k,1,:]) + a2


#  Hydrogen for Chair graphene
for j in range(Nsuper1):
    for k in range(Nsuper2):
        rh[j,k,0,:] = rh1 + (j-1)*a1new + (k-1)*a2new
        rh[j,k,1,:] = rh2 + (j-1)*a1new + (k-1)*a2new
        rh[j,k,2,:] = (rh[j,k,0,:]) + a2
        rh[j,k,3,:] = (rh[j,k,1,:]) + a2



#print (#no atom type numbers for POSCAR...just a list of positions
for j in range(Nsuper1):
    for k in range(Nsuper2):
        for iat in range(Ncarbon):
            file1.write('%8.6f %8.6f %8.6f \n' % ( r[j,k,iat,0], r[j,k,iat,1], r[j,k,iat,2]))

    #print hydrogen atoms
    for j in range(Nsuper1):
        for k in range(Nsuper2):
            for iat in range(NH):
#                 if ~((j == 1 && k==1 && iat ==1)||(j == 1 && k==1 && iat ==3)...
#                         ||(j == 1 && k==2 && iat ==1) ||(j == 1 && k==2 && iat ==3))
                file1.write('%8.6f %8.6f %8.6f \n' % ( rh[j,k,iat,0], rh[j,k,iat,1], rh[j,k,iat,2]))
#                 end


# Adatoms
# rAd =(r(1,1,1,:))' + [0 0 dAd]
rAd = [0, 0, dAd] # Hollow site
file1.write('%8.6f %8.6f %8.6f \n' % ( rAd[0], rAd[1], rAd[2]))
# rAd =(r(1,1,3,:))' + [0 0 dAd]
# file1.write('%8.6f %8.6f %8.6f \n' % ( rAd(1), rAd(2), rAd(3))            
# rAd =(r(1,2,1,:))' + [0 0 dAd]
# file1.write('%8.6f %8.6f %8.6f \n' % ( rAd(1), rAd(2), rAd(3))      
# rAd =(r(1,2,3,:))' + [0 0 dAd]
# file1.write('#8.6f #8.6f #8.6f \n' % ( rAd(1), rAd(2), rAd(3))

file1.close()