#! /usr/bin/env python

#Computes polarization state, and the DN states are obtained by Dmol, but with an arbitrary
#q in each shell.  So six parameters are Vo, ro, Vo, ro, qPol, qDmol
#Put same vo and ro in both ground and excited shells
#CAUTION:  !!!!!!!!!!!!!!!!!!!!!!!!!!!!  q in dmol is the charge to have in the states in the end, after the +2 ion is assumed for generating the states!

############## Initialize ##############

# Starting directory is optimize/loops for now


#print sys.path
import sys, commands, os, time
sys.path.append( '/fslhome/apps/lib/numpy-1.0.4/lib64/python2.3/site-packages/numpy/core' )
sys.path.append( '/fslapps/lib/numpy-1.0.4/lib64/python2.3/site-packages/' )
sys.path.append( '/fslapps/lib/scipy-0.6.0/lib64/python2.3/site-packages' )

#from  os import *
from numpy import *
from  minfunc_s2p2ddmol import *

# Printing options
set_printoptions(precision=3, suppress=True, linewidth=120)
set_printoptions(threshold=nan)

# minfunctrial for testing
def minfunctrial(x): 
 	energy = 5.0*(x[0]-10)**2 + 3.0*(x[1]-2)**2 +2.0*(x[2]-30)**4.00*(x[3]-4)**2+6.0*(x[4]-50)**2+7.0*(x[5] - 6)**2 + 8.0*(x[6] - 7)**2
# 	energy = energy/10000
	return energy

#####################################################
print ' '
print 'Starting optimization loop from masterdmol.py'
print ' '
#####################################################

############## initialize ##############
initial0()
		
############## Hyd initialize ##############			
parametersH = array([30.0, 3.8, 0.0, 0])  #These stay fixed in this version

############## C intialize ##############
parameters = array([100.,  3.5,  250., 4.0, -1, -1, 0]) 
############## Procedure ##############
spacemin = array([100.,  3.5,  250., 4.0, -4.,-4.,  -4])
spacemax = array( [100.,  3.5,  250. , 4.0,  4, 4 ,  4])
str2 = 'Seven parameters are Vo, ro, Vo, ro, qpertDmol, qpertDmol, qpertPol,  (relative to same occupation as unexcited)'
#spacemin = array([-10,-10, -10,-10, -10, -10])
#spacemax = array( [60, 60,60, 60, 60, 60])
summaryfile = open('summary.dat', 'w')
summaryfile2 = open('summary2.dat', 'w')
summaryfile3 = open('summary3.dat', 'w')
summaryfile4 = open('summary4.dat', 'w')
summaryfile.writelines('Benzene iexc 1'+'\n')
summaryfile2.writelines('Benzene iexc 2'+'\n')
summaryfile3.writelines('Diamond iexc 1'+'\n')
summaryfile4.writelines('Diamond iexc 2'+'\n')
str1 = 'iwalk, istep, n,  Vo, ro, Vo, ro, qpertDmol, qpertDmol, qpertPol,  Tfinal, etot, etrans, xcomp , fit'
summaryfile.writelines(str1+'\n')
summaryfile2.writelines(str1+'\n')
summaryfile3.writelines(str1+'\n')
summaryfile4.writelines(str1+'\n')

summaryfile.writelines(str2+'\n')
summaryfile2.writelines(str2+'\n')
summaryfile3.writelines(str2+'\n')
summaryfile4.writelines(str2+'\n')

summaryfile.flush()
summaryfile2.flush()
summaryfile3.flush()
summaryfile4.flush()

sys.stdout.flush()
Nparam = len(parameters)
Nwalks = 500
Nsteps = 3
report = zeros([Nwalks,9], float)
xguess = zeros(Nparam,float)
xlist = zeros([Nsteps,Nparam], float)
energylist = zeros(Nsteps, float)
spread = 0.1
coolfactor = 0*0.005

for j in range (0,Nparam):
	if (parameters[j] > spacemax[j]):
		print 'parameter ' + str(j+1) + ' larger than max specified'
		sys.exit(0)
	if (parameters[j] < spacemin[j]):
		print 'parameter  ' + str(j+1) + ' smaller than min specified'
		sys.exit(0)

# Initialize, make, etc

# Loop over Nwalks 
n=0   # number of function calls
for iwalk in range( 0 , Nwalks ):  # really goes to Nwalks -1....python's boundaries. 
    print 'iwalk', iwalk    
    # Loop over Nstep random search points about guess point
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
                        	trialwidth = spacemax[iparam] - spacemin[iparam] 
				intervalmin = max(spacemin[iparam] , xguess[iparam]  - spread * trialwidth )
                        	intervalmax = min(spacemax[iparam] , xguess[iparam]  + spread * trialwidth )        
                        	xguess[iparam] = (intervalmax + intervalmin)/2.0 + (random.random() - 0.5) * (intervalmax - intervalmin)
#				print 'xguess[iparam]  - spread * trialwidth'
#				print xguess[iparam]  - spread * trialwidth
#				print 'intervalmin, intervalmax' 
#				print intervalmin, intervalmax
#				print 'spread * trialwidth,  (intervalmax + intervalmin)/2.0', spread * trialwidth,  (intervalmax + intervalmin)/2.0
#				print ' (random.random() - 0.5) * (intervalmax - intervalmin)',  (random.random() - 0.5) * (intervalmax - intervalmin)
#

        xlist[istep,:] = xguess
	print 'xlist', xlist
#        energylist[istep] = minfunc6param(xguess)

	functiontest = 0  #if =1, then use numerical function to test parameter distributions
	if ( functiontest == 1 ):
		energylist[istep] = minfunctrial(xguess)
		#print 'function', energylist[istep]

	else:  # not a test
		sys.stdout.flush()
		results = minfunc(parametersH, xguess,n) # call to evaluate create, fireball 
		print 'results', results
		energylist[istep] = results[0][2]
		saveout = sys.stdout                                     
		sys.stdout = summaryfile 
		i1 = 0 
		print  (3*"%4i") % (iwalk, istep, n),  (12*"%7.4f ") % (xguess[0],xguess[1],xguess[2],xguess[3],xguess[4],xguess[5], xguess[6],
		 results[i1][1], results[i1][2], results[i1][3], results[i1][4], results[i1][0]	) #only first benzene result is used to minimize [0]
		sys.stdout.flush()	
		
		sys.stdout = summaryfile2 
		i1 = 1 
		print  (3*"%4i") % (iwalk, istep, n),  (12*"%7.4f ") % (xguess[0],xguess[1],xguess[2],xguess[3],xguess[4],xguess[5], xguess[6],
		 results[i1][1], results[i1][2], results[i1][3], results[i1][4], results[i1][0]	)  
		sys.stdout.flush()	

		sys.stdout = summaryfile3 
		i1 = 2 
		print  (3*"%4i") % (iwalk, istep, n),  (12*"%7.4f ") % (xguess[0],xguess[1],xguess[2],xguess[3],xguess[4],xguess[5], xguess[6],
		 results[i1][1], results[i1][2], results[i1][3], results[i1][4], results[i1][0]	)  
		sys.stdout.flush()	

		sys.stdout = summaryfile4 
		i1 = 3 
		print  (3*"%4i") % (iwalk, istep, n),  (12*"%7.4f ") % (xguess[0],xguess[1],xguess[2],xguess[3],xguess[4],xguess[5], xguess[6],
		 results[i1][1], results[i1][2], results[i1][3], results[i1][4], results[i1][0]	)  
		sys.stdout.flush()	
				
		sys.stdout = saveout  

    print 'xlist'
    print  xlist	
    print 'energylist'
    print energylist
    emin = min(energylist)
    imin = argmin(energylist)
    print 'imin', imin
    xguess = xlist[imin,:]
   
    report[iwalk,0] = iwalk
    report[iwalk,1:8] = xguess[0:7]
    report[iwalk,8] = emin

	

    spread = spread * (1-coolfactor)

# end loop over Nwalks

print 'report'
print report
summaryfile.close()




