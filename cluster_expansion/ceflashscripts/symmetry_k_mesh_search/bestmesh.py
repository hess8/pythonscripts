import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection
    
from unconstrainedSVmin import unconstrainedSVsearch

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,nonzero,float64
fprec=float64
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from itertools import combinations

def bestmeshEigen(Blatt,Nmesh):
    '''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer 3x3 matrix
    So K = B Inv(M) .  Work in the inverse space of this problem, where we can work with M instead of Inv(M). 
    T(InvK) =  T(InvB)T(M).  
    
    Define S = T(InvK), and the real lattice A = T(InvB). So S = A T(M) is a superlattice on the real lattice.
           
    Minimization scheme
    
    Nothing calls this routine'''
    
    ##############################################################
    ########################## Script ############################
   
    M = zeros((3,3),dtype = int)
    S = zeros((3,3),dtype = fprec)
    B = lattice()
    A = lattice()
    K = lattice()
       
    B.vecs = Blatt/2/pi  #Don't use 2pi constants in RL here.
    #############End BCT lattice
    eps = 1.0e-6
    B.det = det(B.vecs)
    B.Nmesh = Nmesh
    print 'B vectors';print B.vecs #
    #print 'B transpose'; print transpose(B.vecs)
    print 'Det of B', B.det
    print 'Orth Defect of B', orthdef(B.vecs)
    print 'Surf/vol of B', surfvol(B.vecs)
    
    [B.symops,B.nops] = getGroup(B.vecs)
    print 'Number of symmetry operations', B.nops
    eigendirs = zeros([3,3,B.nops],dtype = int)
    #print 'symmetry operations of B\n'
    #for j in range(nopsB):
    #    print j
    #    op = array(symopsB[:,:,j])
    #    print op
    #find real lattice
    A.vecs = transpose(inv(B.vecs))
    A.det = det(A.vecs)
    A.Nmesh = Nmesh
    print;print 'A vectors';print A.vecs
    print 'Det of A', A.det
    print 'Orth Defect of A', orthdef(A.vecs)
    print 'Surf/vol of A', surfvol(A.vecs)
    
    [A.symops,A.nops] = getGroup(A.vecs)
    if A.nops != B.nops: 
        sys.exit('Number of operations different for A and B; stop')
    testvecs = [];testindices = []
#    print 'symmetry operations R of A\n'
    for k in range(A.nops):
        print 
        print k
        op = array(A.symops[:,:,k])
        print'symop R of A'; print trimSmall(op)
        m = trimSmall(dot(dot(inv(A.vecs), A.symops[:,:,k]),A.vecs))          
        print'symop m'; print m  
#        print 'det(m)', det(m)              
        'Take eigenvectors in cartesian space'
        [vals,vecs]=eig(op) 
        print 'eigen of m',vals
        print 'eigenvecs are calculated in cartesian space'; print vecs
        #transform to m space
        for i in range(3): vecs[:,i] = dot(inv(A.vecs),vecs[:,i])
        print 'eigenvecs in m space'; print vecs           
        print 'scaled to integers'
        for i in range(3): vecs[:,i] = vecs[:,i]/abs(vecs[:,i])[nonzero(vecs[:,i])].min()
        vecs = rint(vecs)        
        print vecs
        eigendirs[:,:,k]= vecs       
        #find operations with nondegenerate real eigenvalues
        print 'nonDegen', nonDegen(vals)
        for i in nonDegen(vals):
            if not matchDirection(transpose(vecs[:,i]),testvecs): #keep only unique directions    
                testvecs.append(vecs[:,i].real/abs(vecs[:,i])[nonzero(vecs[:,i])].min())
                testindices.append([k,i])
    #print; print oplist;
    print 'testvecs'; print testvecs
    #print testindices
    MT = zeros((3,3),dtype = fprec)
    
    if len(testvecs) == 0:
        print 'No eigen directions'
        [M,K.vecs] = unconstrainedSVsearch(B)
        if det(K.vecs)==0:
            sys.exit('Det(K) is zero after unconstrained search! Stop')
        if not checksymmetry(K.vecs,B):
            sys.exit('Symmetry missing in mesh! Stop')
    #    MT = unconstrainedmin(B.vecs)
    if len(testvecs) == 1:
        print 'Only 1 eigen direction'
        #Choose this one and the other two in the plane perpendicular to this. 
        MT[:,0] = testvecs[0]
#       print 'testindices',testindices
        kop = testindices[0][0] #useful operator 
        ieigen = testindices[0][1] #index of the only eigendirection 
        op = array(A.symops[:,:,kop])
    #    print trimSmall(op)
    
#        find one other direction in the plane perp to the eigendireation; either degenerate eigenvalue will do.
        otherindices = nonzero(array([0,1,2])-ieigen)
        print eigendirs[:,:,otherindices[0][0]]
        MT[:,1] = eigendirs[:,:,kop][:,otherindices[0][0]]
        #Make 3rd vector perp as possible to the other two 
        ur0 = dot(A.vecs,MT[:,0])/norm(dot(A.vecs,MT[:,0])) #unit vectors in real space
        ur1 = dot(A.vecs,MT[:,1])/norm(dot(A.vecs,MT[:,1]))
        ur2 = cross(ur0,ur1)
        print ur0
        print ur1
        print ur2
        print 'ur2 transformed to m space'; print dot(inv(A.vecs),ur2)
        mvec = dot(inv(A.vecs),ur2) #transformed to M space, but real
        mvec = rint(mvec/abs(mvec[nonzero(mvec)]).min()) #Scale so smallest comp is 1
        MT[:,2] = mvec       
        print 'MT from single operator';print MT
        print 'starting superlattice'; print dot(A.vecs,MT)
        
    #    Q2 = MT2mesh_three_ns(MT,B)
        Q2 = MT2mesh(MT,B,A)
        if checksymmetry(Q2,B):
            SV = surfvol(Q2)
    #        print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'  
            K.vecs = Q2                
        else:
            print'Q from single operator fails symmetry'    
    
    if len(testvecs) == 2:
        print 'Only 2 eigen directions'
        MT[:,0] = testvecs[0]
        MT[:,1] = testvecs[1]
        #Make 3rd vector perp as possible to the other two 
        ur0 = dot(A.vecs,MT[:,0])/norm(dot(A.vecs,MT[:,0])) #unit vector in real space
        ur1 = dot(A.vecs,MT[:,1])/norm(dot(A.vecs,MT[:,1]))
        ur2 = cross(ur0,ur1)
        MT[:,2] = rint(dot(inv(A.vecs),ur2))
        print 'MT from two eigen directions';print MT
    #    Q2 = MT2mesh_three_ns(MT,B)
        Q2 = MT2mesh(MT,B)
        if checksymmetry(Q2,B):
            SV = surfvol(Q2)
            print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'  
            K.vecs = Q2                
        else:
            print'Q fails symmetry'  
                        
    if len(testvecs) >= 3:
        print 'MT from three eigen directions'
        testvecstrials = [list(x) for x in combinations(testvecs,3)]
        print testvecstrials    
        bestindex = -1 
        bestcost = 1000 
        for i,vecs in enumerate(testvecstrials):
            print; print 'trial',i
            print vecs
            MT[:,0] = vecs[0]
            MT[:,1] = vecs[1]
            MT[:,2] = vecs[2]
            print 'MT'; print MT
            print 'det MT', det(MT)
            if not areEqual(det(MT),0):
                Q2 = MT2mesh(MT,B)
                if checksymmetry(Q2,B):
                    Nscale =1*.8; Ncost = Nscale * abs((B.det/det(Q2))-B.Nmesh)/B.Nmesh 
                    cost = surfvol(Q2)*(1+Ncost)
                    print cost
                    if cost<bestcost: 
                        bestcost = cost; 
                        bestindex = i; 
                        K.vecs = Q2
                    print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'                  
                else:
                    print'Q from trial %i fails symmetry' % i
        print '___________ Best mesh ___________'
        print 'trial', bestindex
    if checksymmetry(K.vecs,B):
        print K.vecs
        K.det = abs(det(K.vecs))
        print 'N of mesh', B.det/K.det
        SV = surfvol(K.vecs)
        print round(surfvol(K.vecs),4),round(orthdef(K.vecs),4),'SV of Q2,','OD' 
    else:
        print'K mesh fails symmetry'    