#! /usr/bin/env python

############## Initialize ##############

# Starting directory is optimize/loops for now


#print sys.path
import sys, commands, os, time
sys.path.append( '/fslhome/apps/lib/numpy-1.0.4/lib64/python2.3/site-packages/numpy/core' )
sys.path.append( '/fslapps/lib/numpy-1.0.4/lib64/python2.3/site-packages/' )
sys.path.append( '/fslapps/lib/scipy-0.6.0/lib64/python2.3/site-packages' )

#from  os import *
from numpy import *
from  minfunc_s2p2d import *

# minfunctrial for testing
set_printoptions(precision=3, suppress=True, linewidth=120)
def minfunctrial(x): 
 	energy = 5.0*x[0]**2 + 3.0*x[1]**4 + 7.0*x[0]*x[1]**2 +2.0*x[2]+5.0*x[3]+6.0*x[4]+7.0*x[5]
# 	energy = energy/10000
	return energy

#####################################################
print ' '
print 'Starting optimization loop from masterloop.py'
print ' '
#####################################################

############## initialize ##############
initial0()
		
############## Hyd initialize ##############			
parametersH = array([30.0, 3.8, 0.0, 0])  #These stay fixed in this version

############## C intialize ##############
parameters = array([0,  3.018, 0,  4.185, 300., 0.5]) 
############## Procedure ##############
summaryfile = open('summary.dat', 'w')
str1 = 'iwalk, istep, n,  parameters 1-6,  Tfinal, etot, etrans, xcomp , fit'
summaryfile.writelines(str(str1)+'\n')
summaryfile.flush()
sys.stdout.flush()
Nparam = len(parameters)
Nwalks = 500
Nsteps = 3
report = zeros([Nwalks,8], float)
xguess = zeros(Nparam,float)
xlist = zeros([Nsteps,Nparam], float)
energylist = zeros(Nsteps, float)
spacemin = array([0,1.0, 0, 1.0, 0, 1.0])
spacemax = array( [300, 8.0, 300, 8.0, 300, 10.0])
spread = 0.2
coolfactor = 0*0.005


# Initialize, make, etc

# Loop over Nwalks 
n=0   # number of function calls
for iwalk in range( 0 , Nwalks ):  # really goes to Nwalks -1....python's boundaries. 
    print 'iwalk', iwalk    
    # Loop over Nstep random search points about guess point
    print 'range', range( 0, Nsteps )
    for istep in range( 0, Nsteps ):
	n = n + 1
	print 'istep', istep         
        # make new random poin
	if (iwalk == 0) and (istep == 0):
	    #first point, which is n = 0
                 xguess = parameters
	else:
		for iparam in range( 0 , Nparam ):
			#if iparam % 2 <> 0: #remainder: only the odd (starting from zero) ones are changed, to leave Vo alone
                        	intervalmin = max(spacemin[iparam] , xguess[iparam] * (1.0 - spread ))
                        	intervalmax = min(spacemax[iparam] , xguess[iparam] * (1.0 + spread ))        
                        	xguess[iparam] = (intervalmax + intervalmin)/2.0 + (random.random() - 0.5) * (intervalmax - intervalmin)
                    	#end
                #end
        #end
#	print xguessx(
#	print xlist
        xlist[istep,:] = xguess
	print 'xlist', xlist
#        energylist[istep] = minfunc6param(xguess)

#	print 'function', minfunctrial(xguess)
	sys.stdout.flush()
	results = minfunc(parametersH, xguess,n)
	print 'results', results
	energylist[istep] = results[0]


#	energylist[istep] = minfunctrial(xguess)
#	print 'energylist', energylist


#	Write results to file for each function call (not just chosen min)

	#	summaryline = [iwalk, istep, xguess[0:6], results[1:5], results[0]]



	saveout = sys.stdout                                     
	sys.stdout = summaryfile  
	print  (3*"%3i") % (iwalk, istep, n),  (11*"%6.3f ") % (xguess[0],xguess[1],xguess[2],xguess[3],xguess[4],xguess[5], 
		 results[1], results[2], results[3], results[4], results[0]	) 
	sys.stdout.flush()	
	sys.stdout = saveout  

#	summaryfile.writelines(str(summaryline)+'\n')

        # end loop over Nstep
    # find point with lowest energy, make new guess point
    emin = min(energylist)
    imin = argmin(energylist)
    xguess = xlist[imin,:]
    report[iwalk,0] = iwalk
    report[iwalk,1:7] = xguess[0:6]
    report[iwalk,7] = emin

	

    spread = spread * (1-coolfactor)

# end loop over Nwalks

print 'report'
print report
summaryfile.close()



