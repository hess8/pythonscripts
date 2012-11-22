def initial0(): 

    """ Initializes, makes, etc"""

    ##############  intialize ##############30            
    import os
#    from  os import *
    run1=os.popen('./initial0_s2p2d.bash', 'r')
#    os.wait()
    print 'Done initial0'
    return run1.read()

    
#    return ''

def initial1(parameters,n): 

    """ Initializes atom 1 info needed to find minimum energy (or other value), called from masterloop"""

    ############## Hyd intialize ##############30            
    import os
    from numpy import *
#    print 'Initial 1 started'
    infile = '../beginrun/initial.run.inputHDN'
    outfile = '../beginrun/initial.run.input'
    in1 = open(infile, 'r')
    out1 = open(outfile, 'w')
    lines = in1.readlines()
    for j in range(0,4): #Vo, ro go in lines 11- (python 10-)
        lines[10+j] = str( parameters[j])  + '\n'
    for j in range(0,len(lines)): 
        out1.write(str(lines[j]))
    in1.close()
    out1.close()
    run1 = os.popen('./initial1_s2p2d.bash ' + str(n), 'r')
    #os.wait()
    print 'Done initial1'
    return run1.read()

def initial2(parameters,n): 

    """ Initializes atom 2 info needed to find minimum energy (or other value), called from masterloop"""
    ############## C intialize ##############30
    import os
    from numpy import *
    infile = '../beginrun/initial.run.input.sp.splitPol'
    outfile = '../beginrun/initial.run.input'
    in1 = open(infile, 'r')
    out1 = open(outfile, 'w')
    lines = in1.readlines()
    lines[15] = str( parameters[6])  + '\n'  # 7th parameter element, qpert, goes in line 16
    #lines 7 and 8 contain normal q.  qpert Dmol are added to these for atoms 1,2. They go in lines 17,18
    lines[16] = str(float(lines[6]) + parameters[4]) + '\n'
    lines[17] = str(float(lines[7]) + parameters[5]) + '\n'
    for j in range(0,4): #Vo, ro go in 22-25 (for dmol choice)
        lines[21 + j] = str( parameters[j])  + '\n'
    
    for j in range(0,4): #Vo, ro go in 26-29 excited states
        lines[25 + j] = str( parameters[j])  + '\n'
    for j in range(0,len(lines)): 
        out1.write(str(lines[j])) # write out all lines, including changed ones. 
    
        
    in1.close()
    out1.close()

    run1 = os.popen('./initial2_s2p2d.bash ' + str(n), 'r')
    #os.wait()
    print 'Done initial2'
    return run1.read()

def minfunc(parameters1, parameters2, n): 

    """ Prepares atomic files, runs Create and FB and reads items of interest"""
    import os, sys
    from numpy import *
#    import initial1
#    import initial2
    print 'Evaluating minimization function'
    # Prepare H
    print 'Initializing 1'
    initial1(parameters1,n)
    print 'Initializing 2'
    # Prepare C
    initial2(parameters2,n)
    print 'Run create and FB'
    print 'n = ', n
    sys.stdout.flush()
    # Run create and FB
    run1 = os.popen('./runCreateFB.bash ' + str(n), 'r')
    os.wait()
    print 'Create and FB done'

    # Get data from result: benzene
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
    results1 = [fit, Tfinal, etot, etrans, xcomp]      #Only benzene1 results are used for minimization
        
    # Get data from result: benzene2
    file = '../fireballrun2/fb'
    Tcheck=os.popen('awk -f awk_Tfinal ' + file , 'r')
    Tfinal=float(Tcheck.read())
    Tcheck.close()
    print 'T-final2', Tfinal

    etotcheck=os.popen('awk -f awk_etot ' + file , 'r')
    etot=float(etotcheck.read())
    print 'etot2', etot
    etotcheck.close()

    file = '../fireballrun2/Eparts.dat'
    etranscheck=os.popen(' awk -f awk_etrans ' + file , 'r')
    etrans=float(etranscheck.read())
    print 'etrans2', etrans
    etranscheck.close()

    file = '../fireballrun2/answer.bas'
    xcompcheck=os.popen('awk -f awk_answer ' + file , 'r')
    xcomp=float(xcompcheck.read())
    xcompcheck.close()
    print 'xcomp2', xcomp
    
    ideal=1.2091
    diff = abs(xcomp - ideal)
    fit = diff/1e-2 +  (etot + 86.1)/0.1 + Tfinal/5
    
    results2 = [fit, Tfinal, etot, etrans, xcomp]

        # Get data from result: diamond1
    
    file = '../diamond1/fb'
    Tcheck=os.popen('awk -f awk_Tfinal ' + file , 'r')
    Tfinal=float(Tcheck.read())
    Tcheck.close()
    print 'T-final D1', Tfinal

    etotcheck=os.popen('awk -f awk_etot ' + file , 'r')
    etot=float(etotcheck.read())
    print 'etot D1', etot
    etotcheck.close()

    etranscheck=os.popen(' awk -f awk_gap ' + file , 'r')
    etrans=float(etranscheck.read())
    print 'etrans D1', etrans
    etranscheck.close()

    file = '../diamond1/answer.bas'
    xcompcheck=os.popen('awk -f awk_answer ' + file , 'r')
    xcomp=float(xcompcheck.read())
    xcompcheck.close()
    print 'xcomp D1', xcomp

    ideal = -0.44625000  # diamond
    diff = abs(xcomp - ideal)
    fit = diff/3e-3 +  (etot + 157.8)/0.1 + Tfinal/5 - 6.4
    
    resultsDmnd1 = [fit, Tfinal, etot, etrans, xcomp]


        # Get data from result: diamond2
    
    file = '../diamond2/fb'
    
    Tcheck=os.popen('awk -f awk_Tfinal ' + file , 'r')
    Tfinal=float(Tcheck.read())
    Tcheck.close()
    print 'T-final D2', Tfinal

    etotcheck=os.popen('awk -f awk_etot ' + file , 'r')
    etot=float(etotcheck.read())
    print 'etot D2', etot
    etotcheck.close()

    etranscheck=os.popen(' awk -f awk_gap ' + file , 'r')
    etrans=float(etranscheck.read())
    print 'etrans D2', etrans
    etranscheck.close()

    file = '../diamond2/answer.bas'
    xcompcheck=os.popen('awk -f awk_answer ' + file , 'r')
    xcomp=float(xcompcheck.read())
    xcompcheck.close()
    print 'xcomp D2', xcomp

    ideal = -0.44625000  # diamond
    diff = abs(xcomp - ideal)
    fit = diff/3e-3 +  (etot + 157.8)/0.1 + Tfinal/5
    
    resultsDmnd2 = [fit, Tfinal, etot, etrans, xcomp]

    results = [results1, results2,resultsDmnd1, resultsDmnd2 ]


    return results

