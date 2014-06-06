'''Convention here is COLUMNS of matrices as vectors'''
################# functions #######################
from numpy import array, cos, sin,arccos, dot, cross, pi,  floor, sum, sqrt, exp, log, asarray
from numpy import sign, matrix, transpose,rint,inner,multiply,size,argmin,argmax,round,ceil
from numpy import zeros,nonzero,float64, sort, argsort, mod, amin, amax
fprec=float64
#from numpy.matlib import zeros, matrix #creates np.matrix rather than array, but limited to 2-D!!!!  uses *, but array uses matrixmultiply
from numpy.linalg import norm, det, inv, eig
from numpy import int as np_int
from random import random, randrange
from ctypes import byref, cdll, c_double, c_int
import time, os, subprocess, sys

utilslib =  cdll.LoadLibrary('/fslhome/bch/vaspfiles/src/hesslib/hesslib.so') 
#had to copy and rename Gus's routine to the one below because ctypes could never find the one with the right name
getLatticePointGroup = utilslib.symmetry_module_mp_get_pointgroup_

def readfile(filepath):
    file1 = open(filepath,'r')
    lines = file1.readlines()
    file1.close()
    return lines

def writefile(lines,filepath):
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

def getline(index,filepath): 
    lines = readfile(filepath)
    return lines[index].strip()
    

def readstructs(filepath):
    lines = readfile(filepath)
    structs = [line.split()[1] for line in lines] #2nd column
    return structs

def readlatt(filepath):
    lines = readfile(filepath)
    latt = zeros((3,3),dtype=float)
    latt[:,0] = array(lines[2].split()[0:3])
    latt[:,1] = array(lines[3].split()[0:3])
    latt[:,2] = array(lines[4].split()[0:3])
    return latt

def natoms_zeros(filepath): #if there is a zero in the natoms list, remove it
    lines = readfile(filepath)
    os.system('rm POSCAR')
#    print 'n atoms line', lines[5]
    atoms = lines[5].split()
    if '0' in atoms: atoms.remove('0')
    lines[5] = ' '.join(atoms)+'\n'
#    print 'after', lines[5]
    writefile(lines,filepath)
    
def scaleposcar(scale2):
    '''replaces 'scale factor' with the numeric value'''
    lines = readfile('POSCAR')
#    os.system('rm POSCAR')
    lines[1]=lines[1].replace('scale factor', str(scale2))
    writefile(lines,'POSCAR')
    print 'POSCAR scale', scale2
    
def getscale(atomic, structchar):
    commstr= 'aconvasp --proto=%s:%s | aconvasp --poscar>POSCAR' % (structchar+'1',atomic)   #for POSCAR creation of first in list, to get scale      
    os.system(commstr)
    file1 = open('POSCAR','r')
    poscar = file1.readlines()
#    poscar = nstrip(poscar)
    scale = poscar[1].strip().split()[0] #string
    return float(scale)

def getscalePOSCAR(): #assumes a POSCAR with volume factor of 1 is already in the main directory, with the correct scale
    file1 = open('POSCAR','r')
    poscar = file1.readlines()
#    poscar = nstrip(poscar)
    scale = poscar[1].strip().split()[0] #string
    return float(scale)

def getL(platt):
    A = readlatt('POSCAR')
    print 'Lattice'; print A
    L = dot(inv(platt),A)
    return L 

def writekpts_cubic_n(n,shift):
    file1 = open('KPOINTS','w')
    file1.write('BCH equiv kpts'+'\n') 
    file1.write('0\n')
    file1.write('Cartesian\n')
    file1.write('%20.15f 0.00 0.00\n' % float(1.0/n))
    file1.write('0.00 %20.15f 0.00\n' % float(1.0/n))
    file1.write('0.00 0.00 %20.15f\n' % float(1.0/n))
    file1.write('%s %s %s\n' % (shift, shift, shift))
    file1.close()

def writekpts_fcc_n(n,shift):
    file1 = open('KPOINTS','w')
    file1.write('BCH equiv kpts'+'\n') 
    file1.write('0\n')
    file1.write('Cartesian\n')
    file1.write('0.00 %20.15f %20.15f\n' % (0.5*float(1.0/n),0.5*float(1.0/n)))
    file1.write('%20.15f 0.00 %20.15f\n' % (0.5*float(1.0/n),0.5*float(1.0/n)))
    file1.write('%20.15f %20.15f 0.00 \n' % (0.5*float(1.0/n),0.5*float(1.0/n)))
    file1.write('%s %s %s\n' % (shift, shift, shift))
    file1.close()  
    
def makestr2poscar(struct):
        print os.system('makestr.x %s %s' % ('../../struct_enum.out', struct)) #creates poscar-like vasp.0xxxx file
        makestr_out = subprocess.check_output('ls vasp.0*', shell=True).strip()
        os.system('cp %s POSCAR' % makestr_out)
        natoms_zeros('POSCAR')  
    