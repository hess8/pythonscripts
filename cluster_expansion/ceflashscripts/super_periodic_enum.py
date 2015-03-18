#!/usr/bin/python
'''
Analyzes two struct_enum.out files:  one that has the superperiodic structures removed, 
and one that has kept them. Finds the extras structures in the 2nd, and renumbers them
to the end of a copy of the first.
The ordering in the files seem scrambled, so have to do a line-by-line search
'''

import sys,os,subprocess
from numpy import zeros,transpose,array,sum,float64,rint,mean
from ceroutines import readfile,writefile

maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1/enum_superperiodic'
os.chdir(maindir)
os.system('cp ../enum/struct_enum.out struct_enum_part.out')
os.system('cp struct_enum_part.out struct_enum.out')
linesfull = readfile('{}/struct_enum_full.out'.format(maindir))
linespart = readfile('{}/struct_enum_part.out'.format(maindir))
linesout = []

outfile = open('{}/struct_enum.out'.format(maindir),'a')

superlines = []

startCol = 61 #drop first columns (chars) that have N and degeneracy, count at particular size, etc info.  

for i,lfull in enumerate(linesfull[17:]):
    print i
    restFull = lfull[startCol:]  
    match = False
    for ipart,lpart in enumerate(linespart[17:]):  
        restPart = lpart[startCol:] 
        if restFull == restPart: #this is not superperiodic
            match = True
#            print 'match'
            break           
    if not match: #superperiodic
        superlines.append(lfull[15:])
        linespart
#Write out superperiodic structures
print 'len', len(superlines)
nPartMax = int(linespart[-1].split()[0])  #from last line...maximum N
for i,line in enumerate(superlines):
    outfile.write('{:12d}   {}'.format(nPartMax+1+i,line))     
outfile.close    
print 'Done'

