import os,subprocess
maindir = '/bluehome/bch/vasprun/graphene.structures/'
newdir = maindir +'transmet.half_graphane/'
#dirs = ['half_graphane/','h.half_graphane2.1/']
#typedirs = ['initial_relax/','final_relax/']

dirs = ['/fslhome/bch/vasprun/graphene.structures/transmet.half_graphane/h.half_graphane2.1/dos/']

elements = [
'Co','Cr','Fe','Hf','Ir'
,'Mn','Mo','Nb_pv','Ni','Os',
'Pd','Pt','Re','Rh','Ru','Ta','Tc','Ti','V','W','Zr'
]

os.chdir(maindir)
for dir in dirs:
    for element in elements:
        path = dir + 'adatom_'+element+'/ISPIN_2/'
        print element, path
#        print(['cp',frompath, topath])
        subprocess.call(['rm','-r', path])
    
print "done"   

##for isolated folders:
#frompath2 = '/bluehome/bch/vasprun/graphene.structures/half_graphane/isolated/'
#topath2 = newdir + 'isolated/electron_relax/'
#for element in elements:
#    frompath = frompath2 + 'adatom_'+element+'/'
#    print element, topath2
##        print(['cp',frompath, topath])
#    subprocess.call(['cp','-r', frompath, topath2])
#    
#print "done"   