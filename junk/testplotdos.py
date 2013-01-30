import sys, os #, vasputil_dosplot_fun
from vasputil_dosplot_fun import vasputil_dosplot 
mainDir ='/fslhome/bch/vasprun/graphene.structures/transmet.half_graphane/h.half_graphane2.1/dos/'
atomDir = 'adatom_Co/ISIF_4/IBRION_2/ISPIN_2/'
subdir = 'dos/'
dir = mainDir + atomDir +subdir
vasputil_dosplot([], ["DOSCAR"], dir) #options, args, dir


