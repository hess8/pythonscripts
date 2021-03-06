def plotArray(x,y,matrix1,plotfile1,title1,xlabel1,ylabel1,plotmax):
    '''plots colored matrix for 2 d array'''
#    from __future__ import division
    from matplotlib.patches import Patch
    from pylab import *
    print plotfile1
#    print x, y
#    x=np.append(x,x[-1])#duplicate last value in extra slot, so plot will show all rows/columns 
#    y=np.append(y,y[-1])
#    print x,y  
    X,Y = meshgrid(x, y)
    Z = matrix1
    fig = figure()
    pcolor(X, Y, Z, cmap=cm.hot, vmax = plotmax)
    xlim((x.min(),x.max()))
    title(title1)
    xlabel(xlabel1)
    ylabel(ylabel1)
    colorbar()
    show()
    fig.savefig(plotfile1)

def plotxy(x,y,plotfile1,title1,xlabel1,ylabel1):
#    from __future__ import division 
    from matplotlib.pyplot import *
    fig = figure()
    plot(x, y, 'ro')
    title(title1)
    xlabel(xlabel1)
    ylabel(ylabel1)
#    ylim((-11.5,-10.5))
    show() 
    fig.savefig(plotfile1)  

def collate_plots(checkedList,plotName):
    import sys, os, subprocess
    import numpy as num
    sys.path.append('/fslhome/bch/pythonscripts/pythonscripts_gr/half_graphane_scripts/analysis_scripts')
    from analysisTools import getElement
    
    '''Creates an HTML page with the plots and labels'''
    if not os.path.exists('plots'):
        subprocess.call(['mkdir','plots'])
    nRow = 7  # put in periodic table order, with duplicates at end
    plotType = plotName.split('.')[-2]
    print plotType
    collatefile  = open('%splots.htm' % plotType,'w')
    collatefile.write(' <html>\n <HEAD>\n<TITLE> %s </TITLE>\n</HEAD>\n' % plotType)
    collatefile.write(' <BODY>\n <p style="font-size:20px"> <table border="1">\n <tr>\n') #start row
    iImage = 0        
    for path in checkedList:
        element = getElement('adatom_',path)
        elementPlotName = element + plotName
        try: 
            subprocess.call(['cp','%s%s' % (path,plotName),'plots/'+elementPlotName])     
        except:
            print 'copy plot failed', path      
        if not element in ['Cr_pv', 'Mn_pv']: # don't want these in first table
            iImage += 1            
            collatefile.write('<td><p><img src="plots/%s" ></p><p>%s</p></td>\n' % (elementPlotName,element))#Image and element under it
            if num.mod(iImage,nRow) == 0: 
                collatefile.write('</tr>\n<tr>\n') #end of row, begin new
    collatefile.write(' </tr></table> \n') #end of row and table                
    for path in checkedList:
        #get element name
        element = getElement('adatom_',path)
        elementPlotName = element + plotName              
        if element  in ['Cr_pv', 'Mn_pv']: # put these last            
            collatefile.write('<td><p><img src="plots/%s" ></p><p>%s</p></td>\n' % (elementPlotName,element))#Image and element under it  
    collatefile.write(' </BODY> </html>') #end of file 
    collatefile.close() 


    
def vasputil_dosplot(options, args, dir):
#!/usr/bin/python
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
# Copyright (c) 2008, 2009, 2010 Janne Blomqvist

# This source code file is subject to the terms of the MIT (Expat)
# License. See the file LICENSE for details.

    """Example plotting script demonstrating how to use the vasputil.dos module
    and the ase.calculators.vasp.VaspDos class."""
    import sys, os
    print dir
    os.chdir(dir)
    
    sys.path.append('/fslhome/bch/pythonscripts/ase-3.6.0.2515/')
    sys.path.append('/fslhome/bch/pythonscripts/vasputil-5.7/')
    sys.path.append('/fslhome/bch/pythonscripts/vasputil-5.7/scripts/')
    sys.path.append('/fslhome/bch/pythonscripts/vasputil-5.7/scripts/vasputil')
    
    import vasputil.dos as d
    import ase.calculators.vasp as v
    import matplotlib.pyplot as plt
    from optparse import OptionParser
        
    usage = """%prog [options] DOSCAR atom
 
    atom is the index of the atom, starting from 0. If omitted, plots the total
    DOS, or the integrated total DOS if additionally the -i option is present.
    """
    parser = OptionParser(usage)
    parser.add_option("-f", "--fermi", dest="fermi", help="Fermi energy")
    parser.add_option("-o", "--orbital", dest="orb", help="Orbital, either \
    string or integer. For non-spin polarized s=0, p=1, \
    d=2, if spin polarized s+=0, s-=1 etc. If phase factors: s, py, pz, px, \
    dxy, dyz, dz2, dxz, dx2, and equivalent as previous for spin polarized etc.")
    parser.add_option('-i', '--integrated', dest='integ', action='store_true',
                      help='Show the integrated total DOS')
    (options, argsdummy) = parser.parse_args()

# By default VaspDos constructor tries to read a file called "DOSCAR".
    dc = v.VaspDos(args[0])
    if options.fermi:
        fermi = float(options.fermi)
        dc.efermi = fermi
    if options.orb:
        try:
            orb = int(options.orb)
        except ValueError:
            orb = options.orb
    else:
        orb = 0
    
    en = dc.energy # This is the x-axis in a typical DOS plot.
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if len(args) == 2:
        atom = int(args[1])
        label = "Atom " + str(atom) + " orbital " + str(orb)
        ax.plot(en, dc.site_dos(atom, orb), label=label)
    elif options.integ:
        label = 'Integrated total DOS'
        if len(dc.integrated_dos.shape) == 1:
            ax.plot(en, dc.integrated_dos, label=label)
        else:
            ul = label + ' spin up'
            ax.plot(en, dc.integrated_dos[0, :], label=ul)
            dl = label + ' spin down'
            ax.plot(en, dc.integrated_dos[1, :], label=dl)
    else:
        label = 'Total DOS'
        if len(dc.dos.shape) == 1:
            ax.plot(en, dc.dos, label=label)
        else:
            ul = label + ' spin up'
            ax.plot(en, dc.dos[0, :], label=ul)
            dl = label + ' spin down'
            ax.plot(en, dc.dos[1, :], label=dl)
    d.set_labels()
    plt.show()
    plt.savefig('dos')
    
    # For a more involved plot, create your own plot script or run interactively
    # via 'ipython -pylab'
    # Below are some example pylab commands for plotting
    
    # For good plots in PS or PDF formats begin your script with
    # from matplotlib import rc
    # rc('ps', usedistiller='xpdf')
    # rc('text', usetex=True)
    
    #subplot(211)
    # Get the s and p DOS of atom 1.
    #plot(en, dc.site_dos(1, 0), "k-", label="Al s")
    #plot(en, dc.site_dos(1, 1), "k-.", label="Al p")
    #legend()
    #xlim(-15,5)
    #xticks(arange(0))
    #subplot(212)
    # s-DOS of atom 2. 
    #plot(en, dc.site_dos(2, 0), "k-", label="Al s")
    #xlim(-15,5)
    #ylim(0,.29)
    
    # And to save the figs in EPS/PDF, use
    # d.set_labels()
    # savefig('tex_demo.eps')
    # Fonts in legend touch upper border, test with your version of matplotlib
    # savefig('tex_demo_xpdf.pdf')
    # os.system('epstopdf tex_demo.eps')
    
    # Finally, show picture
    # show()
    
    
########### Colorline function ########## similar to cline in Matlab.    
    '''
Defines a function colorline that draws a (multi-)colored 2D line with coordinates x and y.
The color is taken from optional data in z, and creates a LineCollection.

z can be:
- empty, in which case a default coloring will be used based on the position along the input arrays
- a single number, for a uniform color [this can also be accomplished with the usual plt.plot]
- an array of the length of at least the same length as x, to color according to this data
- an array of a smaller length, in which case the colors are repeated along the curve

The function colorline returns the LineCollection created, which can be modified afterwards.

From http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb 

:

# Example 1: Sine wave colored by time

x = np.linspace(0, 4.*np.pi, 1000)
y = np.sin(x)

fig, axes = plt.subplots()

colorline(x, y) %note z is optional

plt.xlim(x.min(), x.max())
plt.ylim(-1.0, 1.0)
plt.show()



# Example 2 Shifting a sine wave: Example modified from mpltools by Tony S Yu

x = np.linspace(0, 4.*np.pi, 1000)
y = np.sin(x)

fig, axes = plt.subplots()

N = 10
for i in xrange(N):
    color = i / float(N)
    shift = 0.2 * i
    colorline(x-shift, y, color, cmap="cool") %Here color is a single number, so each line is a single color
    
plt.xlim(x.min(), x.max())
plt.ylim(-1.0, 1.0)
plt.show()  



See also: plt.streamplot
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


# Data manipulation:

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments


# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    
    ax = plt.gca()
    ax.add_collection(lc)
    
    return lc
        
    
def clear_frame(ax=None): 
    # Taken from a post by Tony S Yu
    if ax is None: 
        ax = plt.gca() 
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False) 
    for spine in ax.spines.itervalues(): 
        spine.set_visible(False) 

########### end Colorline function ########## 