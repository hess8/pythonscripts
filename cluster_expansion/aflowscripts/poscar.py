from numpy import dot, cross, arccos, zeros, pi
from numpy.linalg import norm
import operator
from ctypes import byref, cdll, c_double

class POSCAR(object):
    """Reads and writes POSCAR files in VASP format."""
    def __init__(self, path = None, lines = None):
        self.name = ""
        self.scale = 0.0
        #A list of the lattice vectors in the POSCAR
        self.avecs = []
        #A list containing the number of each species of atom that will be
        #declared with dvectors after the "D" symbol
        self.types = []
        #True if the coordinates are specified direct (vs. cartesian)
        self.direct = True
        #A list of the direct vector to atom occupation sites in the cell.
        self.atoms = []
        #A list of the concentrations of atomic species for each of the
        #atomic sites listed in self.atoms.
        self.concentrations = []
        #Stoichiometry constraints for each of the atomic sites listed
        #in self.atoms. Since multiple sites are required to have stoichiometry
        #constraints, the same value is repeated for all such sites.
        self.stoichiometry = []

        #Private vars for optimization
        self._dconcentrations = []

        # Length of lattice vectors
        self._lengths = []
        # Angles betwenn lattice vectors
        self._angles = []
        # Reciprocal lattice vectors
        self._bvecs = []
        # Lengths of the reciprocal lattice vectors
        self._bvecs_lengths = []
        

        if path is not None:
            self.read(path)
        if lines is not None:
            self.fromlines(lines)
    
    @property
    def orthogonality_defect(self):
        """Orthogonality defect of real space (avecs)"""
        volume_avecs = abs(dot(cross(self.avecs[0],self.avecs[1]),self.avecs[2]))
        prod_avecs = reduce(operator.mul,self.abc,1)
        volume_bvecs = abs(dot(cross(self.bvecs[0],self.bvecs[1]),self.bvecs[2]))
        prod_bvecs = reduce(operator.mul,self.abc_bvecs,1)
        return prod_avecs/volume_avecs

    @property
    def recip_orthogonality_defect(self):
        """Orthogonality defect of reciprocal space vectors (bvecs)"""
        volume_bvecs = abs(dot(cross(self.bvecs[0],self.bvecs[1]),self.bvecs[2]))
        prod_bvecs = reduce(operator.mul,self.abc_bvecs,1)
        return prod_bvecs/volume_bvecs

    @property
    def od(self):
        """Shorter name to orthogonality defect"""
        return self.orthogonality_defect

    @property
    def rod(self):
        """Shorter name to recip. orthogonality defect"""
        return self.recip_orthogonality_defect

    @property
    def abc(self):
        """Lengths of the 3 lattice vectors"""
        if len(self._lengths) != 3:
            self._lengths = [norm(A) for A in self.avecs]
        return self._lengths

    @property
    def abc_bvecs(self):
        """Lengths of the 3 reciprocal vectors"""
        if len(self._bvecs_lengths) != 3:
            self._bvecs_lengths = [norm(B) for B in self.bvecs]
        return self._bvecs_lengths

    @property
    def lengths(self):
        """Another handle to the lengths of the 3 lattice vectors"""
        return self.abc

    @property
    def angles(self):
        """The 3 angles between the lattice vectors"""
        # Can we set the __str___ so that we can pretty print these?
        if len(self._angles) != 3:
            self._angles = zeros(3)
            self._angles[0] = arccos(dot(self.avecs[1],self.avecs[2])/self.abc[1]/self.abc[2])
            self._angles[1] = arccos(dot(self.avecs[2],self.avecs[0])/self.abc[2]/self.abc[0])
            self._angles[2] = arccos(dot(self.avecs[0],self.avecs[1])/self.abc[0]/self.abc[1])
            self._angles *= 360./2/pi # Output the results in degrees rather than radians
        return self._angles

    @property
    def abg(self):
        """Another handle to the 3 angles between the lattice vectors"""
        return self.angles

    @property
    def rank(self):
        """Gets the number of unique species in the POSCAR definition."""
        return len(self.types)

    @property
    def numsites(self):
        """Returns the number of sites defined in the poscar."""
        return len(self.atoms)

    @property
    def dconcentrations(self, index):
        """Returns a list of concentrations in decimal format such that
        each value represents the probability of atom occupying that site."""
        if type(self._dconcentrations[index]) == type([]):
            return self._dconcentrations[index]            
        else:
            return 1
    @property
    def volume(self):
        return  abs(dot(cross(self.avecs[0],self.avecs[1]),self.avecs[2]))

    @property
    def bvecs(self):
        if len(self._bvecs) != 3:
            self._bvecs = [zeros(3),zeros(3),zeros(3)]
#            print self.volume
#            print "cross: ",cross(self.avecs[1],self.avecs[2])
            self._bvecs[0] = 2*pi/self.volume*cross(self.avecs[1],self.avecs[2])
            self._bvecs[1] = 2*pi/self.volume*cross(self.avecs[2],self.avecs[0])
            self._bvecs[2] = 2*pi/self.volume*cross(self.avecs[0],self.avecs[1])
        return self._bvecs

    def read(self, path):
        """Reads a POSCAR from the specified path and assigns properties."""
        with open(path) as f:
            lines = f.readlines()
        self.fromlines(lines)
            
    def fromlines(self, lines):
        """Assigns POSCAR properties from a list of strings representing
        the standard form of a poscar file."""
        #Get all the information that will never change
        self.name = lines[0]
        self.scale = float(lines[1])
        self.avecs = [[ float(i) for i in j.split() ] for j in lines[2:5]]
        if self.scale < 0:
            volume_avecs = abs(dot(cross(self.avecs[0],self.avecs[1]),self.avecs[2]))
            self.scale = (-self.scale/volume_avecs)**(1/3.0)    
        self.avecs = [[ self.scale*self.avecs[i][j] for j in [0,1,2] ] for i in [0,1,2]]
        self.scale = 1.0
        self.types = [ int(i) for i in lines[5].split() ]
        self.direct = lines[6][0] == "D"
        #Read in the direct vectors for each of the atom types in self.types
        dindex = 7 #The index of the d-vector currently being read in
        #(offset by 7).
        for sindex in range(self.rank):
            scount = self.types[sindex]
            for i in range(scount):
                v = [ float(j) for j in lines[dindex].split()[:3] ]
                a = AtomicSite(v, sindex)
                self.atoms.append(a)
                dindex += 1
                
            

        #Now we can read any additional information about occupation probabilities
        #for atomic sites that are not pure.
        if dindex < len(lines) - 1:
            print 'dindex', dindex
            print  'len(lines)', len(lines)
            for sindex in range(self.rank):
                scount = self.types[sindex]
                for i in range(scount):
                    cline = lines[dindex].split()
                    if cline[0] != "P":
                        #This is not a pure site, read in the concentrations and
                        #optionally the stoichiometry
                        self._parse_concentration(cline)
                    else:
                        #This is a pure site, just store it as such
                        self.concentrations.append("P")
                        self._dconcentrations.append(1)
                        self.stoichiometry.append(None)
                    dindex += 1
        else:
            #They must want pure on all the sites as specified in the types
            for sindex in range(self.rank):
                scount = self.types[sindex]
                for i in range(scount):
                    #This is a pure site, just store it as such
                    self.concentrations.append("P")
                    self._dconcentrations.append(1)
                    self.stoichiometry.append(None)
    
    def write_poscar(self,filename):
        '''BCH.  Doesn't support concentrations '''
        el = '\n'
        file1 = open(filename,'w')
        file1.write(self.name+el)
        file1.write(str(self.scale)+el)
        for i in [0,1,2]:
            file1.write('%20.14f %20.14f %20.14f \n' % (self.avecs[i][0],self.avecs[i][1],self.avecs[i][2]))     
        str1 = ''
        for i in range(len(self.types)):
            str1 = str1 + str(self.types[i])+' '
        file1.write(str1+el)
        if self.direct:
            file1.write('Direct'+el)            
        else:
            file1.write('Error'+el)
        
        str2 = ''
        for i in range(len(self.atoms)):
            file1.write('%20.14f %20.14f %20.14f \n' % (self.atoms[i].vector[0],self.atoms[i].vector[1],self.atoms[i].vector[2]))
#            file1.write( '  '.join(map(str,self.atoms[i].vector)) + el)
        file1.close()       

    def mink_reduce(self,eps):
        """Reduce the basis to the most orthogonal set.
           A Minkowski-reduced (via a "greedy algorithm basis) """
#        utilslib =  cdll.LoadLibrary('/Users/hart/codes/celib/trunk/libutils.so')
        utilslib =  cdll.LoadLibrary('/fslhome/bch/cluster_expansion/theuncle/celib/trunk/libutils.so')
        mink = utilslib.vector_matrix_utilities_mp_minkowski_reduce_basis_
        
        a = load_ctypes_3x3_double(self.avecs)
        b = ((c_double * 3) *3)()
        epsIN = c_double(eps)
        
        mink(byref(a),byref(b),byref(epsIN))
        bout = unload_ctypes_3x3_double(b)
        return bout

    def kmink_reduce(self,eps):
        """Reduce the reciprocal vectors to the most orthogonal set.
           A Minkowski-reduced (via a "greedy algorithm basis) """
        utilslib =  cdll.LoadLibrary('/fslhome/bch/cluster_expansion/theuncle/celib/trunk/libutils.so')
        mink = utilslib.vector_matrix_utilities_mp_minkowski_reduce_basis_
        
        a = load_ctypes_3x3_double(self.bvecs)
        b = ((c_double * 3) *3)()
        epsIN = c_double(eps)
        
        mink(byref(a),byref(b),byref(epsIN))
        bout = unload_ctypes_3x3_double(b)
        return bout

def load_ctypes_3x3_double(IN):
    """Make a 3x3 array into the right thing for ctypes"""
    a = ((c_double * 3) *3)()
    for i in range(3):
        for j in range(3):
            a[i][j] = c_double(IN[i][j])
    return a

def unload_ctypes_3x3_double(OUT):
    """Take a ctypes array and load it into a 3x3 python list"""
    a = zeros((3,3))
    for i in range(3):
        for j in range(3):
            a[i][j] = OUT[i][j]
    return a
    
 
#    def _parse_concentration(self, cline):
#        """Parses out the concentration rules for the occupation at the site
#        specified in the concentration line of the POSCAR file.
#
#         - cline: a list of concentration rules. The first rule deals with
#           probability of occupation. The second (optionally) enforces
#           stoichiometric constraints.
#
#        Concentrations are specified as 1:2:3 implying a 1/6 chance to have
#        species 1 at the site and so on.
#        """
#        #The concentration exists for stoich. cases and normal ones
#        self.concentrations.append([ int(n) for n in cline[0].split(":") ])        
#        #Calculate the decimal probabilities of the concentrations
#        #These are used for random dice selections
#        total = sum(self.concentrations[-1])
#        for n in self.concentrations[-1]:
#            self._dconcentrations.append(sum(self._dconcentrations) + float(n) / total)
#
#        if len(cline) > 1:
#            #We have stoichiometry to worry about
#            self.stoichiometry.append([ int(n) for n in cline[1].split(":") ])
#
class AtomicSite(object):
    """Represents a single atom in a crystal."""
    def __init__(self, vector, species):
        self.vector = vector
        self.species = species
        self.direct = True
        