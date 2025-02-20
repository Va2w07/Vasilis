import pylab
import numpy
import random
import scipy.special
import scipy.stats
import pickle
import copy 

epsilon0 = 8.8e-12
m0 = 9.10938291e-31 #Electron mass
hbar = 6.58211928e-16 #in eV.s
kb = 8.6173324e-5 #in eV/K
me = 0.067*m0 
mh = 0.5*m0
Eg = 1.424 #eV
Temp = 300
Ephoton = 1.55 #eV 800nm
Ephonon = 0.035
w0 = Ephonon/hbar
q = 1.60217657e-19
fb = 1./(numpy.e**(Ephonon/(kb*Temp))-1)
e0 = 12.95*epsilon0
einf = 10.9*epsilon0
epsilonp = 1./(1./einf-1./e0)
M = -me*mh/(me-mh)

class Carrier:
    def __init__(self,i):
        self.charge = -1*q
        self.x = 0
        self.y = 0
        self.z = 0
        self.kx = 0
        self.ky = 0
        self.kz = 0
        self.t0 = 0
        self.E = 0
        self.tf = 0
        self.tfold = 0
        self.destroy = False
        self.i = i
        self.m = me
        
    def drift(self, t):
        self.kx += self.E*self.charge*t/(hbar*q)
        return
    
    def scatter(self):
        rGamma = gamma*numpy.random.random()
                
        sr = recombinationRate()/gammascale
        if rGamma <= (sr):
            self.destroy = True
            return
        
        if self.charge < 0:
            spe =  phononRate(self.kx, self.m, 1)/gammascale

            if numpy.isnan(spe+sr):
                spe = 0
            if rGamma <= (spe):
                self.kx = phononNewK(self.kx, self.m, 1)
                return

            spa =  phononRate(self.kx, self.m, -1)/gammascale
            if numpy.isnan(spa):
                spa = 0
            if rGamma <= (spe+spa+sr):
                self.kx = phononNewK(self.kx, self.m, -1)
                return
    
            si = impurity_scatter(self.kx, self.m)/gammascale
            if numpy.isnan(si):
                si = 0
            if rGamma <= (spe+spa+sr+si):
                self.kx = impNew(self.kx)
                return            
        return
        
def Ek(k,m):
    return hbar**2*k**2*q/(2*m)

ni = 1e15 #Defect Density

def d_scatter(k,m, n = 1e19):
    # n will need to be changed once we have messing to calculate the local carrier densities
    k = abs(k)
    kd = numpy.sqrt(q*n/(epsilon0*kb*Temp))
    return 2*ni*q**2*k/(numpy.pi*hbar**3*epsilon0*kd**2*(kd**2+4*k**2))

def kmax(k, m, aore):
    return k*abs(1+(numpy.sqrt(1-aore*Ephonon/Ek(k,m))))

def kmin(k, m, aore):
    return k*abs(1-(numpy.sqrt(1-aore*Ephonon/Ek(k,m))))

def phononRate(k, m, aore): #aore = -1 for absorption and +1 for emission
    k = abs(k)
    return m*w0/(4*numpy.pi*epsilonp*1*k*hbar**2)*(fb+0.5+aore*0.5)*numpy.log(kmax(k,m,aore)/kmin(k,m,aore))

def recombinationRate():
    return 1./(200e-15)


dscat = []
pescat = []
pascat = []
rescat = []

for E in numpy.arange(0, 0.70, 0.00001):
    k = numpy.sqrt(2*me*E*q)/(hbar*q)
    
    dscat.append(d_scatter(k,me))
    pescat.append(phononRate(k,me,1))
    pascat.append(phononRate(k,me,-1))
    rescat.append(recombinationRate())

pylab.plot(numpy.arange(0, 0.70, 0.00001), dscat)
pylab.plot(numpy.arange(0, 0.70, 0.00001), pescat)
pylab.plot(numpy.arange(0, 0.70, 0.00001), pascat)
pylab.plot(numpy.arange(0, 0.70, 0.00001), rescat)

pylab.yscale('log')
pylab.savefig('scatterRates2.pdf', format='pdf')

pylab.show()
1/0
carrierList = []
with open( "carFinal", "rb" ) as file:
    while 1:
        try:
            carrierList.append(pickle.load(file))
        except EOFError:
            break


carrierList[0] = sorted(carrierList[0], key=lambda car: car.kx) #Put in generation time order
Kbine = numpy.arange(-8.7e8, 8.8e8, 1e7)*0
Kbinh = numpy.arange(-8.7e8, 8.8e8, 1e7)*0
Xbine = numpy.arange(-170e-6, 170e-6, 1e-6)*0
Xbinh = numpy.arange(-170e-6, 170e-6, 1e-6)*0
Ebine = numpy.arange(0, 0.70, 0.01)*0
Klist = []
ktot=0
Ee = []


for carr in carrierList[0]:
    kx = 87+carr.kx/1e7
    x = 170+carr.x*1e6
    E = Ek(carr.kx, carr.m)*100
    if kx>87*2:
        continue
    if carr.charge < 0:
        Klist.append(kx)
        Kbine[int(kx)] += 1
        ktot += abs(kx-87)
        Xbine[int(x)] += 1
        Ee.append(Ek(carr.kx, carr.m))
        Ebine[int(E)]+=1
    else:
        Kbinh[int(kx)] += 1
        Xbinh[int(x)] += 1
print ktot/sum(Kbine)

pylab.bar(numpy.arange(0, 0.70, 0.01), Ebine, width = 0.01)
pylab.savefig('finalDistEg.pdf', format='pdf')
pylab.show()
pylab.clf()

pylab.bar(numpy.arange(-170e-6, 170e-6, 1e-6),  numpy.array(Xbine), width = 1e-6, color = 'b')
pylab.bar(numpy.arange(-170e-6, 170e-6, 1e-6),  numpy.array(Xbinh), width = 1e-6, color = 'r', alpha=0.5)
pylab.savefig('finalXg.pdf', format='pdf')
pylab.show()
pylab.clf()

pylab.bar(numpy.arange(-8.7, 8.8, 0.1),Kbine, width = 0.1)
pylab.bar(numpy.arange(-8.7, 8.8, 0.1),Kbinh, width = 0.1, color = 'r', alpha = 0.5)
pylab.xlim(-6,6)
pylab.savefig('finalDistkxg.pdf', format='pdf')
pylab.show()
