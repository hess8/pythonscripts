from numpy import cos,sin,pi

#a1 =    array([2.13129, -1.2305, 0 ])     
#a2 =    array([2.13129, 1.2305, 0  ])    
#a3s =    array([.00000000,  0.000000000, 15.00]) 
##this one will be unchanged in supercell
#theta = arccos(dot(a1,a2)/norm(a1)/norm(a2)) #in case of distorted hexagonal lattice
#thetadeg = theta*180/pi
## Doubling size of unit cell to allow boat structures. New primitive cell can be made of a diagonal and a horizontal vector
## a1new = a1 + a2#5# horizontal...need to double number of atoms in unit cell.
## a2new = a2 - a1 #(gives a vertical lattice vector)
#a1new = a1 # Don't change
#a2new = a2 # 
##  Atom positions from vasp primitive cell relaxation
#rc1 =      array([1.46545,      0.00000,      0.22856])        
#rc2 =      array([2.93090,      0.00000,     -0.22856])         
#rh1 =      array([1.46545,      0.00000,      1.33850])      
#rh2 =      array([2.93090,      0.00000,     -1.33850 ])   


print cos(pi/3)
print sin(pi/3)
print 1/2.0/cos(pi/6)
print 1/2.0/cos(pi/6)*2
print 1/2.0/cos(pi/6)*3
print 2.13129/1.2305