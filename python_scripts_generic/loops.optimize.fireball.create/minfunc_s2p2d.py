def initial0(): 

	""" Initializes, makes, etc"""

	##############  intialize ##############30			
	import os
#	from  os import *
	run1=os.popen('./initial0_s2p2d.bash', 'r')
	os.wait()
	print 'Done initial0'
	return run1.read()

	
#	return ''

def initial1(parameters,n): 

	""" Initializes atom 1 info needed to find minimum energy (or other value), called from masterloop"""

	############## Hyd intialize ##############30			
	import os
	from numpy import *
	print 'Initial 1 started'
	paramfile = open('../beginrun/paramatom', 'w')
	print 'opened paramatom 1', 
	for j in range(0,len(parameters)):
		paramfile.writelines(str(parameters[j])+'\n')
	paramfile.close()		
	run1 = os.popen('./initial1_s2p2d.bash ' + str(n), 'r')
	os.wait()
	print 'Done initial1'
	return run1.read()

def initial2(parameters,n): 

	""" Initializes atom 2 info needed to find minimum energy (or other value), called from masterloop"""
	############## C intialize ##############30
	import os
	from numpy import *
	paramfile = open('../beginrun/paramatom', 'w')
	for j in range(0,len(parameters)):
		paramfile.writelines(str(parameters[j])+'\n')
	paramfile.close()
	run1 = os.popen('./initial2_s2p2d.bash ' + str(n), 'r')
	os.wait()
	print 'Done initial2'
	return run1.read()

def minfunc(parameters1, parameters2, n): 

	""" Prepares atomic files, runs Create and FB and reads items of interest"""
	import os
	from numpy import *
#	import initial1
#	import initial2
	print 'Evaluating minimization function'
	# Prepare H
	print 'Initializing 1'
	initial1(parameters1,n)
	print 'Initializing 2'
	# Prepare C
	initial2(parameters2,n)
	print 'Run create and FB'
	print 'n = ', n
	# Run create and FB
	run1 = os.popen('./runCreateFB.bash ' + str(n), 'r')
	os.wait()
	print 'Create and FB done'

	# Get data from result
	Tcheck=os.popen('awk -f awk_Tfinal ../fireballrun/fb', 'r')
	Tfinal=float(Tcheck.read())
	Tcheck.close()
	print 'T-final', Tfinal

	etotcheck=os.popen('awk -f awk_etot ../fireballrun/fb', 'r')
	etot=float(etotcheck.read())
	print 'etot', etot
	etotcheck.close()

	etranscheck=os.popen(' awk -f awk_etrans ../fireballrun/Eparts.dat', 'r')
	etrans=float(etranscheck.read())
	print 'etrans', etrans
	etranscheck.close()

	xcompcheck=os.popen('awk -f awk_answer ../fireballrun/answer.bas', 'r')
	xcomp=float(xcompcheck.read())
	xcompcheck.close()
	print 'xcomp', xcomp

	ideal=1.2091
	diff = abs(xcomp - ideal)
	fit = diff/1e-2 +  (etot + 86.1)/0.1 + Tfinal/5
	
	results = [fit, Tfinal, etot, etrans, xcomp]

	return results

