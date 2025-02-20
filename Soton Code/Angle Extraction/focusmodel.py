import numpy
import pylab

L = 20e-3
n1 = 1
n2 = 1.43

theta = numpy.arange(0,0.3,0.01)
theta2 = numpy.arcsin(numpy.sin(theta)*n1/n2)

z = L*(numpy.tan(theta) + numpy.tan(theta2))/numpy.tan(theta)

w0 = 3e8/(1e12*numpy.pi*theta)
w1 = w0*numpy.sqrt(1+(z*numpy.pi*1e12*theta/(3e8))**2)

T1 = 2*n1*numpy.cos(theta)/(n1*numpy.cos(theta2)+n2*numpy.cos(theta))
T2 = 2*n2*numpy.cos(theta2)/(n2*numpy.cos(theta)+n1*numpy.cos(theta2))

T = T1*T2

pylab.figure()
pylab.plot(theta, T/T[0])

pylab.show()
