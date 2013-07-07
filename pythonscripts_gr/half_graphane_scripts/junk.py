import os,subprocess
maindir = '/bluehome/bch/vasprun/graphene.structures/'
newdir = maindir +'transmet.half_graphane/'
dirs = ['h.half_graphane2.1/']
typedirs = ['final_relax/']

elements = ['Co','Cr','Cr_pv','Fe','Hf','Ir','Mn','Mn_pv','Mo','Nb_pv','Ni','Os','Pd','Pt','Re','Rh','Ru','Ta','Tc','Ti','V','W','Zr']
#elements = ['Co']
os.chdir(maindir)
for dir in dirs:
    frompath1 = maindir + dir
    topath1 = newdir + dir
    for type in typedirs:
        frompath2 = frompath1 + type
        topath2 = topath1 + type        
        for element in elements:
            frompath = frompath2 + 'adatom_'+element+'/'
            topath = topath2
            print element, topath2
    #        print(['cp',frompath, topath])
            subprocess.call(['cp','-r', frompath, topath])  
print "done"   
