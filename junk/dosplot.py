#import vasputil_dos
import os, sys

mainDir ='/fslhome/bch/vasprun/graphene.structures/transmet.half_graphane/h.half_graphane2.1/dos/'
atomDir = 'adatom_Co/ISIF_4/IBRION_2/ISPIN_2/'
subdir = 'dos/'
dir = mainDir + atomDir +subdir
print dir
os.chdir(dir)
#vdos = vasputil_dos.LDOS()
#vdos.read_doscar()
#print dos
sys.path.append('/fslhome/bch/pythonscripts/python-ase-3.6.0.2515/')
sys.path.append('/fslhome/bch/pythonscripts/vasputil-5.7/')
sys.path.append('/fslhome/bch/pythonscripts/vasputil-5.7/scripts/')
sys.path.append('/fslhome/bch/pythonscripts/vasputil-5.7/scripts/vasputil')
import vasputil_dosplot