import numpy
import pylab

ang = []
dx = 1e-6 
xlist = numpy.arange(50e-6, 1e-3, dx)
for x in xlist:
    ang.append(numpy.arccos(x/(x+2e-6))*360/(2*numpy.pi))
    
pylab.plot(xlist, ang)
pylab.show()