import pylab
import numpy
import random
import scipy.special
import scipy.stats
import pickle
import copy 

epsilon0 = 8.8e-12
m0 = 9.10938291e-31 #Electron mass
hbar = 1.0545e-34
kb = 1.38e-23
q = 1.60217657e-19

me = 0.067*m0 
mh = 0.5*m0

Eg = 1.424 #eV
Temp = 300

Ephoton = 1.55*q #eV 800nm
Ephonon = 0.035*q

w0 = Ephonon/hbar
fb = 1./(numpy.e**(Ephonon/(kb*Temp))-1)
e0 = 12.95*epsilon0
einf = 10.9*epsilon0
epsilonp = 1./(1./einf-1./e0)
M = -me*mh/(me-mh)

def Ek(k,m):
    return k**2*hbar**2/(2*m)
def kmax(k, m, aore):
    return k*abs(1+(numpy.sqrt(1-aore*Ephonon/Ek(k,m))))

def kmin(k, m, aore):
    return k*abs(1-(numpy.sqrt(1-aore*Ephonon/Ek(k,m))))

def phononRate(k, m, aore): #aore = -1 for absorption and +1 for emission
    k = abs(k)
    return q**2*w0*k/(8*numpy.pi*epsilonp*Ek(k,m))*(fb+0.5+aore*0.5)*numpy.log(kmax(k,m,aore)/kmin(k,m,aore))

dscat = []
pescat = []
pascat = []
rescat = []

for E in numpy.arange(0, 0.70, 0.01):
    E *= q
    k = numpy.sqrt(2*me*E)/(hbar)
    
    pescat.append(phononRate(k,me,1))
    pascat.append(phononRate(k,me,-1))

pylab.plot(numpy.arange(0, 0.70, 0.01), pescat)
pylab.plot(numpy.arange(0, 0.70, 0.01), pascat)

pylab.yscale('log')
pylab.show()