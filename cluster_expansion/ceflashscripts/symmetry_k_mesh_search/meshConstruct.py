'''
Builds a mesh from the symmetry operations and random starting points:
1. nCoarse random points in the cell parallelpiped.  
2. Define cutoff as rc = (nCoarse/Vcell)^1/3
3. Make all copies due to translation and symmetry operations
4. In neighboring cells (3x3x3 construction of plpd), keep only those that 
are within e.g. 1.5 * rc of the cell boundaries
5. If two points are less than e.g. 0.8*rc away from each other, combine them 
into a point at their center. 
5.5 if we have too many delete those that are closest to their neighbors
5.6 If we have too few, get the Delaunay tetrahedra, and add to the center 
of the largest tets. 
6. Keep track of independent points and symmetry partners
7. Calculate the force on each ind. point from all the points that lie within
a distance rForce.  Use a power law: Fi = k*Sum(r_vec_i,j/ rij^p
8.  Move the ind. points in a vector proportional to the force
9.  Copy all the ind. points by trans and pt symmetry again. 
10. repeat from 7 until the force on each point is less than minForce
11. Copy all the ind. points by trans and pt symmetry again. 
________
Fine trahedral mesh
________
1. Find Delaunay tetrahedra
2. Subdivide each tet that has an independent point in it with a formula that 
minimizes the S/V of the 6 new tets. Call these independent children, 
even though they might have parents that are not independent
3. Copy these by trans and point sym, removing duplicates
4. Subdivide the indChildren (repeat 2 and 3) until we have the number of points we want.
5. By symmetry, find the number of truly independent among all points, and their weights. 
6. Write these to KPOINTS: https://cms.mpi.univie.ac.at/vasp/vasp/Entering_all_k_points_explicitly.html 
until we have the number of points we want. 

Disadvantage: nK grows by 6 each time...is that a problem? Not for now. 


'''