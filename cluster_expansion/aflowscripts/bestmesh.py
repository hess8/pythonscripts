import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import svmesh, svmesh1freedir, lattice_vecs, lattice, surfvol, \
    orthdef, icy, isinteger, isequal, isreal, isindependent, trimSmall, cosvecs,  \
    load_ctypes_3x3_double, unload_ctypes_3x3_double, unload_ctypes_3x3xN_double, \
    getGroup, checksymmetry, nonDegen, MT2mesh, matchDirection
    
from nconstrainedSVmin import unconstrainedSVsearch

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import transpose,rint,inner,multiply,size,argmin,nonzero
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from itertools import combinations

def bestmesh(Blatt,Nmesh):
    '''The kmesh can be related to the reciprocal lattice B by  B = KM, where M is an integer 3x3 matrix
    So K = B Inv(M) .  Work in the inverse space of this problem, where we can work with M instead of Inv(M). 
    T(InvK) =  T(InvB)T(M).  
    
    Define S = T(InvK), and the real lattice A = T(InvB). So S = A T(M) is a superlattice on the real lattice.
           
    Minimization scheme'''
    
    ##############################################################
    ########################## Script ############################
    
    #natoms = 3
    #nkppra = 10000
    #nk = int(nkppra/natoms)
    
    #print 'Target N kpoints', Nmesh
    
    M = zeros((3,3),dtype = int)
    S = zeros((3,3),dtype = float)
    B = lattice()
    A = lattice()
    K = lattice()
    
    ##############BCT lattice
    alat = 2*sqrt(5)
    ca = 7/3.
    clat = alat*ca
    #B.vecs = array((  
    #  [   -alat/2,  alat/2,   alat/2],
    #  [   alat/2,  -alat/2,   alat/2],
    #  [   clat/2,   clat/2,   -clat/2]
    #  ), dtype=float)
    
    
    #B.vecs = array((  #C axis along x !####
    #  [   -clat/2,  clat/2,   clat/2],
    #  [   alat/2,  -alat/2,   alat/2],
    #  [   alat/2,   alat/2,   -alat/2]
    #  ), dtype=float)
    
    #print 'B vectors before inverse and transpose';print B.vecs
    #B.vecs = trimSmall(transpose(inv(B.vecs)))
    #############End BCT lattice
    
    ############## Any lattice
    
    #crystal = [1,1,sqrt(2),90,90,120] # [a,b,c,alpha,beta,gamma]
    crystal = [2,3*sqrt(3),11*sqrt(2),90,90,20] # [a,b,c,alpha,beta,gamma]
    #crystal = [2,2,2,80,80,80] # [a,b,c,alpha,beta,gamma]
    
    B.vecs = transpose(lattice_vecs(crystal))
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
    #print 'symmetry operations of B\n'
    #for j in range(nopsB):
    #    print j
    #    op = array(symopsB[:,:,j])
    #    print op
    #find real lattice
    A.vecs = trimSmall(transpose(inv(B.vecs)))
    A.det = det(A.vecs)
    A.Nmesh = Nmesh
    print;print 'A vectors';print A.vecs
    print 'Det of A', A.det
    print 'Orth Defect of A', orthdef(A.vecs)
    print 'Surf/vol of A', surfvol(A.vecs)
    
    [A.symops,A.nops] = getGroup(A.vecs)
    if A.nops != B.nops: 
        sys.exit('Number of operations different for A and B; stop')
    
    print 'symmetry operations R of A\n'
    testvecs = [];testindices = []
    for k in range(A.nops):
        print; print k
        op = array(A.symops[:,:,k])
        print trimSmall(op)
        m = trimSmall(dot(dot(inv(A.vecs[:,:]), A.symops[:,:,k]),A.vecs[:,:])  ) 
        [vals,vecs]=eig(m); vecs = array(vecs)
        print 'symop m'; print m
    
        print 'det(m)', det(m)
        print 'eigen of m',vals
        print vecs
        print 'as rows, scaled to integers'
        print vecs[:,0]/abs(vecs[:,0])[nonzero(vecs[:,0])].min()
        print vecs[:,1]/abs(vecs[:,1])[nonzero(vecs[:,1])].min()
        print vecs[:,2]/abs(vecs[:,2])[nonzero(vecs[:,2])].min()
        
        #find operations with nondegenerate real eigenvalues
    #    print 'nonDegen', nonDegen(vals)
        for i in nonDegen(vals):
            if not matchDirection(transpose(vecs[:,i]),testvecs): #keep only unique directions    
                testvecs.append(vecs[:,i].real/abs(vecs[:,i])[nonzero(vecs[:,i])].min())
                testindices.append([k,i])
    #print; print oplist;
    #print testvecs
    #print testindices
    MT = zeros((3,3),dtype = int)
    
    if len(testvecs) == 0:
        print 'No eigen directions'
        K.vecs = unconstrainedSVsearch(B)
        if det(K.vecs)==0:
            sys.exit('Det(K) is zero after unconstrained search! Stop')
        if not checksymmetry(K.vecs,B):
            sys.exit('Symmetry missing in mesh! Stop')
    #    MT = unconstrainedmin(B.vecs)
    if len(testvecs) == 1:
        print 'Only 1 eigen direction'
        #Choose this one and the other two in the plane perpendicular to this. 
        #Since all symmetry operators will be diagonal in this mesh representaion
        #of eigenvectors of , 
        #
        MT[:,0] = testvecs[0]
        k = testindices[0][0]
        op = array(A.symops[:,:,k])
    #    print trimSmall(op)
        m = trimSmall(dot(dot(inv(A.vecs[:,:]), A.symops[:,:,k]),A.vecs[:,:])  ) 
        [vals,vecs]=eig(m); vecs = array(vecs)
    #    print vecs
        otherindices = nonzero(array([0,1,2])-k)
        v1 = vecs[otherindices[0][0]]
        MT[:,1] = v1/abs(v1)[nonzero(v1)].min()
        #Make 3rd vector perp as possible to the other two 
        ur0 = dot(A.vecs,MT[:,0])/norm(dot(A.vecs,MT[:,0])) #unit vector in real space
        ur1 = dot(A.vecs,MT[:,1])/norm(dot(A.vecs,MT[:,1]))
        ur2 = cross(ur0,ur1)
        MT[:,2] = rint(dot(inv(A.vecs),ur2))
        print 'MT from single operator';print MT
    #    Q2 = MT2mesh_three_ns(MT,B)
        Q2 = MT2mesh(MT,B)
        if checksymmetry(Q2,B):
            SV = surfvol(Q2)
    #        print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'  
            K.vecs = Q2                
        else:
            print'Q fails symmetry'    
    
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
        testvecstrials = [list(x) for x in combinations(testvecs,3)]
    #    print testvecstrials    
        bestindex = -1 
        bestcost = 1000 
        for i,vecs in enumerate(testvecstrials):
            print; print 'trial',i
            MT[:,0] = vecs[0]
            MT[:,1] = vecs[1]
            MT[:,2] = vecs[2]
            print 'MT from three eigen directions';print MT
    #        print 'det MT', det(MT)
            if not isequal(det(MT),0):
                Q2 = MT2mesh(MT,B)
                if checksymmetry(Q2,B):
                    Nscale =1*.8; Ncost = Nscale * abs((B.det/det(Q2))-B.Nmesh)/B.Nmesh 
                    cost = surfvol(Q2)*(1+Ncost)
                    if cost<bestcost: bestcost = cost; bestindex = i; K.vecs = Q2
    #                print round(surfvol(Q2),4),round(orthdef(Q2),4),'SV of Q2,','OD'                  
                else:
                    print'Q fails symmetry'
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