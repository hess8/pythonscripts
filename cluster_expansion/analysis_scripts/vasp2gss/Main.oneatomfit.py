'''
Makes clusters for model fitting. Starts with completed structures.in files, fits and performs ground state search.  

Make changes to lat.in in /enum/
'''

import os, subprocess,sys
from random import seed
from numpy import zeros
from copy import deepcopy


if __name__ == '__main__':
    #Make model changes in enum/lat.in
#    dir = '/fslhome/bch/cluster_expansion/graphene/tm_row1/Cr_pv/fits'
#    os.chdir(dir)
    dir = os.getcwd()   
    os.chdir('../../enum')
    os.system('uncle 10')
    os.chdir(dir)
    os.system('cp ../../enum/lat.in .')
    os.system('rm gss.out')    
    os.system('cp ../../enum/clusters.out .')
    os.system('uncle 15')
    os.system('uncle 21')
    os.system('mv gss.out ../gss/')
    os.chdir('../gss')
    os.system('gnuplot gss_plot.gp')
    
    print 'Done'

    

        
    
    
 
 
 
 
 
 
 
    