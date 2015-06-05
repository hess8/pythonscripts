import os, subprocess, sys, time 

#sys.path.append('/bluehome2/bch/pythonscripts/Erik_scripts/') 
#import plotGraphene

structs = [7755,2446,218,3099,657,1239]


maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1/'
os.chdir(maindir)
dir = '/fslhome/bch/cluster_expansion/graphene/tm_row1/struct_plots'

if not os.path.isdir(dir): os.mkdir(dir)
os.chdir(dir)
os.system('ln -s ~/pythonscripts/Erik_scripts/plotGraphene.py ')

for struct in structs:
    os.system('makestr.x ../enum/struct_enum.out {}'.format(struct))
    
dirList = os.listdir(dir)
for item in dirList:
    for struct in structs:        
        if str(struct) in item:
            print 'plotting structure '+ str(struct)
#            print 'python plotGraphene.py . {}  -u'.format(item)
            os.system('python plotGraphene.py . {}  -u -z'.format(item))
print 'done'