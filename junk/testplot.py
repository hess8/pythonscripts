#import numpy
#import pylab
#import matplotlib
#t = numpy.arange(0.0, 1.0+0.01, 0.01)
#s = numpy.cos(2*2*numpy.pi*t)
#plot(t, s)
#show()
#
#
#
#pylab.plot(t, s)
# 
#pylab.xlabel('time (s)')
#pylab.ylabel('voltage (mV)')
#pylab.title('About as simple as it gets, folks')
#pylab.grid(True)
#pylab.savefig('simple_plot')
#
#pylab.show()
#print matplotlib.__version__

#from numpy import *
#from scipy import *
#from matplotlib import *
#from pylab import *
#
#ion()
#plot([1,2,3],[1,2,3])
#figure(1)
#raw_input(':')
#close(1)

import pylab as p
import numpy as n
p.ion()
p.figure(1)
a = n.array([1,2,3])
p.clf()
p.plot(a)
p.show()


