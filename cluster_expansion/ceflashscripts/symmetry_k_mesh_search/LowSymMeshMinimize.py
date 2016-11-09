'''Routines for minimizing S/V without regard to symmetry.  Useful for lowest symmetry crystals.'''
'''NOT USED FOR NOW...SEE UnconstrainedSVmin.py'''


import os, subprocess, sys, time 

sys.path.append('/bluehome2/bch/pythonscripts/cluster_expansion/aflowscripts/')
from kmeshroutines import lattice,surfvol, orthdef

from kmeshroutines import lattice_vecs, lattice, surfvol, orthdef, icy, isinteger, areEqual, isreal, isindependent, trimSmall, cosvecs

from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import matrix, transpose,rint,inner,multiply,size,argmin,argmax,nonzero,shape
from numpy import zeros #use arrays, not "matrix" class
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy import int as np_int
from numpy import float as np_float
from random import random, randrange, randint
from ctypes import byref, cdll, c_double, c_int

#from numpy import array, arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, matrix, transpose,rint,inner,multiply
#from numpy import zeros as array_zeros
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
#from numpy.linalg import norm, det, inv
#from numpy import int as np_int
#from numpy import float as np_float
#from random import random, randrange
#from ctypes import byref, cdll, c_double, c_int

utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') 
#had to copy and rename Gus's routine to the one below because ctypes could never find the one with the right name

def searchmin(S,A):
    '''MT is transpose(M) '''
    print '----------------------'
    default_length = rint(A.Nmesh**(1/3.0))
    MT = dot(inv(A.vecs),S)
    MTold = MT
    #determine how many full rows (known vectors) there are in S:
    if norm(MTold[:,0])==0: knownvecs = 0; MT[0,0] = default_length # Make first trial vector along "mx"
    elif norm(MTold[:,0])!=0 and norm(MTold[:,1])==0: knownvecs = 1 #Make 2nd trial vector 
    else: knownvecs = 2 
    if norm(MT[:,1])==0: #make arbitrary 2nd trial M vector nearly normal to first   
        closestM_axis = argmax(MT[:,0])#find which axis the M vector is most inclined towards
        #choose arbitrary next axis in cyclic way, and make a unit vector in the plane of MT0 and (next axis)
        nextdir = zeros((3,1),dtype = int); nextdir[icy(closestM_axis,1)] = 1 
        nextdir = default_length*nextdir
        print MT
        print 'nextdir',nextdir 
        print 'norm2', norm(MT[:,0])**2
        orthvec = nextdir - dot(transpose(MT[:,0]),nextdir)[0,0]*MT[:,0]/norm(MT[:,0])**2
        print 'orth',orthvec
        unitvec = orthvec/norm(orthvec)
        print'unit',unitvec
        MT[:,1] = rint(default_length*unitvec)
    if norm(MT[:,2])==0: #make arbitrary 3nd trial M vector nearly normal to first two   
        orthvec = transpose(cross(transpose(MT[:,0]),transpose(MT[:,1])))
        unitvec = orthvec/norm(orthvec)
        MT[:,2] = rint(default_length*unitvec)
#    print 'MT', MT
#    print 'det', det(MT) 
    #find the direction of steepest slope in the cost
    maxsteps = 10000
    istep = 1
    while istep<maxsteps:
        bestindex = changewhich(MT,knownvecs,A)
        print 'bestindex',bestindex            
        if bestindex[1]==0:#found minimum at previous M
            newcost = cost(MT,A)        
            break
        else:
#            print 'value in MT to change', bestindex[0][0], MT[bestindex[0]]
            MT[bestindex[0][0],bestindex[0][1]] += bestindex[1]
            newcost = cost(MT,A)
#            oldcost = newcost
        istep += 1
    #    sys.exit('stop')        
    if istep < maxsteps:
        S = dot(A.vecs,MT)
        K = lattice();K.vecs = inv(transpose(S)); K.det = det(K.vecs)        
        print
        print 'Found minimum after %i steps' % istep
        print 'Symmetry error', round(symmetryErr(S,A),4)       
        print 'An optimum transpose(M):'; print MT
        print 'Number of mesh points', det(S)/A.det
        print 'An optimum superlattice S:'; print S
        print 'An optimum K mesh\n', K.vecs  
        print 'Orth defect',orthdef(K.vecs)
        print 'Surface/vol', surfvol(K.vecs)    
        print 'Mesh vector lengths:'; print norm(K.vecs[:,0]),norm(K.vecs[:,1]),norm(K.vecs[:,2])
#        k0 = K.vecs[:,0]; k1 = K.vecs[:,1]; k2 = K.vecs[:,2]
#        cosgamma = k0.T*k1/norm(k0)/norm(k1)
#        cosalpha = k1.T*k2/norm(k1)/norm(k2)
#        cosbeta =  k2.T*k0/norm(k2)/norm(k0)       
#        print 'Mesh vector cosines:'; print cosalpha, cosbeta, cosgamma
#        print 'Check B = KM   \n', K.vecs*M  
#        print '\n\n\nTranspose for use in VMD or POSCAR'
#        print 'B'; print B.vecs.T
#        print 'K'; print K.vecs.T
        print 
    else:
        print 'Ended without minimum after maximum %i steps' % istep
    return MT


def cost(MT,A):
    if det(MT) == 0: print MT;print sys.exit('Error det(MT)=0; stop in cost()')
    S = dot(A.vecs,MT)
    Nscale =1* .8; Ncost = Nscale * abs((det(S)/A.det -A.Nmesh))/A.Nmesh
    symmScale = 10;  symmCost = symmScale* symmetryErr(S,A)  
#    print 'Ncost, symmCost', Ncost, symmCost
#    print 'for MT:'
#    print MT
    cost = surfvol(S)*(1+Ncost + symmCost)        
    return cost

def symmetryErr(latt,parentlatt):
    '''Returns an error that shows how far a lattice is from obeying the symmetry of the parent lattice'''
    err = 0.0
    for iop in range(parentlatt.nops):
        lmat = array(latt)
        if det(lmat) == 0:print lmat; sys.exit('Determinant is zero; stop')               
        mmat = dot(dot(inv(lmat),parentlatt.symops[:,:,iop]),lmat)
#        print 'mmat in symmetryErr';print mmat
        for i in range(3):
            for j in range(3):
                err += abs(rint(mmat[i,j])-mmat[i,j])
    return err/parentlatt.nops 