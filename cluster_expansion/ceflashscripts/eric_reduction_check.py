import os, subprocess
from ctypes import byref, cdll, c_double
import numpy as np
from numpy import linalg as la
    
def mink_reduce(a,eps):
    """Reduce the basis to the most orthogonal set.
       A Minkowski-reduced (via a "greedy algorithm basis) """
#        utilslib =  cdll.LoadLibrary('/fslhome/eswens13/hess_research/libutils.so')
    utilslib =  cdll.LoadLibrary('/fslhome/eswens13/hess_research/libutils.so')
    ared =((c_double * 3) *3)()
    mink = utilslib.vector_matrix_utilities_mp_minkowski_reduce_basis_ 
    mink(byref(load_ctypes_3x3_double(a)),byref(ared),byref(c_double(eps)))
    ared2 = unload_ctypes_3x3_double(ared)
#    print ared2
#    print ared[1]
    
    return ared2

def load_ctypes_3x3_double(IN):
    """Make a 3x3 array into the right thing for ctypes"""
    a = ((c_double * 3) *3)()
    for i in range(3):
        for j in range(3):
            a[i][j] = c_double(IN[i][j])
    return a

def unload_ctypes_3x3_double(OUT):
    """Take a ctypes array and load it into a 3x3 python list"""
    a = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            a[i][j] = OUT[i][j]
    return a

def orthogonality_defect(a):
        """Orthogonality defect of real space (avecs)"""
        volume_avecs = abs(np.dot(np.cross(a[0],a[1]),a[2]))
        prod_avecs = abc(a)[0]*(abc(a)[1]*abc(a)[2])
       # volume_bvecs = abs(dot(cross(self.bvecs[0],self.bvecs[1]),self.bvecs[2]))
       # prod_bvecs = reduce(operator.mul,self.abc_bvecs,1)
        return prod_avecs/volume_avecs

def abc(a):
        """Lengths of the 3 lattice vectors"""
	b = np.array((la.norm(a[0]), la.norm(a[1]), la.norm(a[2])), dtype = float)
	return b	
		
#a = np.array(([1,.5,1.5], [.5,0,.5],[.5,.5,0]), dtype=float)
#eps = 1e-4
#ared = np.array(mink_reduce(a,eps))

#n = np.dot(ared, a)

#print "a:\n", a
#print "original o_defect: ", orthogonality_defect(a)
#print "\nared:\n", ared
#print "reduced o_defect: ", orthogonality_defect(ared)
#print "\nn\n", n


for i in range(1,51):
	os.chdir('/fslhome/eswens13/fsl_groups/hessgroup/temp/AlIr/f' + str(i))
	
	outfilename = "f" + str(i) + "out"

# Open the original file and write to new one
	try: 
		infile = open('/fslhome/eswens13/fsl_groups/hessgroup/temp/AlIr/f' + str(i) + '/POSCAR', 'r')	
		vec_list = infile.readlines()
	except IOError:
		pass
	infile.close()

	eps = 1e-4
	
	outfile = open('/fslhome/eswens13/hess_research/AlIr/' + outfilename, 'w')

# Write only the lines containing the lattice vectors
	outfile.write(vec_list[2])
	outfile.write(vec_list[3])
	outfile.write(vec_list[4])

	outfile.close()

# Reformat so there is one space between each number
	subprocess.call(['sed','-i', 's/   / /g','/fslhome/eswens13/hess_research/AlIr/'+ outfilename])
	subprocess.call(['sed','-i', 's/  -/ -/g','/fslhome/eswens13/hess_research/AlIr/'+ outfilename])

# Take out the space at the beginning of the line and write to
# final_outs directory	
	try:
		infile2 = open('/fslhome/eswens13/hess_research/AlIr/f' + str(i) + 'out', 'r')
		
	except IOError:
		pass
	
	outfile2 = open('/fslhome/eswens13/hess_research/final_outs/' + outfilename, 'w')
	for line in infile2:
		line = line[1:]
		outfile2.write(line)
	
	infile2.close()	
	outfile2.close()

# Removes all files in the AlIr directory.
	# os.system('rm -rf ~/hess_research/AlIr/*')

# Reads the files into a list of floats and makes workable matrices with them.
	infile3 = open('/fslhome/eswens13/hess_research/final_outs/f' + str(i) + 'out', 'r')
	outfile3 = open('/fslhome/eswens13/hess_research/matrix_out/f' + str(i) + 'out', 'w')
	
	num_list = []
	for line in infile3:
		num_list = num_list + line.split()

	vec1 = [float(num_list[0]), float(num_list[1]), float(num_list[2])]
	vec2 = [float(num_list[3]), float(num_list[4]), float(num_list[5])]
	vec3 = [float(num_list[6]), float(num_list[7]), float(num_list[8])]

	a = np.array(vec1,vec2,vec3)
	ared = np.array(mink_reduce(a, eps))

	outfile3.write('Original: ')
	outfile3.write(str(vec1))
	outfile3.write(str(vec2))
	outfile3.write(str(vec3))

	outfile3.write('\nOriginal Defect:  ' + str(orthogonality_defect(a)))

	outfile3.write('Reduced:  ')
	outfile3.write(str(ared[0]))
	outfile3.write(str(ared[1]))
	outfile3.write(str(ared[2]))
	
	outfile3.write('\nReduced Defect:  ' + str(orthogonality_defect(ared)))	

	infile3.close()
	outfile3.close()