import numpy
import pylab

angrange = numpy.arange(0,0.3,0.01)

n1 = 1
n2 = 1.43
L = 2e-2

dx = L * (1 - (n1 * numpy.cos(angrange) ) / numpy.sqrt(n2**2 - n1**2*numpy.sin(angrange)**2) )

pylab.plot(angrange,dx/dx[0])
pylab.show() 