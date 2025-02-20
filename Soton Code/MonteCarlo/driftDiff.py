import copy
import numpy
import os
import pickle
import scipy.stats

L = 100
N = 1e6

gamma = 1
gammascale = 1e14

T = 0
Tmax = 1e-12
Deltat = 1e-15

sigmax = 1e-6
sigmat = 10e-15

epsilonr = 13.1
epsilon0 = 8.8e-12
m0 = 9.10938291e-31  # Electron mass
hbar = 6.58211928e-16  # in eV.s
kb = 8.6173324e-5  # in eV/K
me = 0.0067 * m0 
mh = 0.5 * m0
Eg = 1.424  # eV
Temp = 300
Ephoton = 1.55  # eV 800nm
Ephonon = 0.035
w0 = Ephonon / hbar
q = 1.60217657e-19
fb = 1. / (numpy.e ** (Ephonon / (kb * Temp)) - 1)
e0 = 12.95 * epsilon0
einf = 10.9 * epsilon0
epsilonp = -e0 * einf / (einf - e0)
M = -me * mh / (me - mh)

xMax = 5e-6
xMin = -5e-6
dx = 10e-9

meshEmpty = numpy.arange(xMin, xMax, dx) * 0
mesh = numpy.arange(xMin, xMax, dx) * 0
pmesh = numpy.arange(xMin, xMax, dx) * 0

ke = 1. / (4 * numpy.pi * epsilon0)

class Carrier:
    def __init__(self, i):
        self.charge = -1 * q
        self.x = 0
        self.y = 0
        self.z = 0
        self.kx = 0
        self.ky = 0
        self.kz = 0
        self.t0 = 0
        self.E = 0
        self.Ex = 0
        self.tf = 0
        self.tfold = 0
        self.destroy = False
        self.i = i
        self.m = me
        
    def drift(self, t):
        self.kx += self.Ex * self.charge * t / (hbar * q)
        return
    
    def scatter(self):
        rGamma = gamma * numpy.random.random()
                
        sr = recombinationRate() / gammascale
        if rGamma <= (sr):
            self.destroy = True
            return
        
        if self.charge != 0:
            spe = phononRate(self.kx, self.m, 1) / gammascale

            if numpy.isnan(spe + sr):
                spe = 0
            if rGamma <= (spe):
                self.kx = phononNewK(self.kx, self.m, 1)
                return

            spa = phononRate(self.kx, self.m, -1) / gammascale
            if numpy.isnan(spa):
                spa = 0
            if rGamma <= (spe + spa + sr):
                self.kx = phononNewK(self.kx, self.m, -1)
                return
    
            si = d_scatter(self.kx, self.m) / gammascale
            if numpy.isnan(si):
                si = 0
            if rGamma <= (spe + spa + sr + si):
                self.kx = impNew(self.kx)
                return            
        return
    
    def freeflight(self):
        numpy.random.seed()
        self.tfold = self.tf
        self.tf = -numpy.log(numpy.random.random()) / gammascale
        return
    
    def propergate(self, t):
        self.x += t * q * hbar * self.kx / self.m + q * self.Ex * t ** 2 / (2 * self.m)
        return

carrierReserve = [Carrier(i) for i in numpy.arange(0, N)]

carrierList = []

DeltaE = Ephoton - Eg
Ei = DeltaE * mh / (me + mh)
 
# Set the maxwell-Boltzmann distribution mean
mean = Ei
maxw = scipy.stats.maxwell(0, 1)
meanp = maxw.stats()[0]
maxw = scipy.stats.maxwell(0, mean / meanp)

ni = 1e15  # Defect Density m^-3
Z = 1

def recombinationRate():
    return 1. / (250e-15)

def kmax(k, m, aore):
    return k * abs(1 + (numpy.sqrt(1 - aore * Ephonon / Ek(k, m))))

def kmin(k, m, aore):
    return k * abs(1 - (numpy.sqrt(1 - aore * Ephonon / Ek(k, m))))

def phononRate(k, m, aore):  # aore = -1 for absorption and +1 for emission
    k = abs(k)
    return m * w0 / (4 * numpy.pi * epsilonp * k * hbar ** 2) * (fb + 0.5 + aore * 0.5) * numpy.log(kmax(k, m, aore) / kmin(k, m, aore))

def phononNewK(k, m, aore):
    if k < 0:
        a = -1
    else:
        a = 1
    return abs(K(Ek(k, m) - aore * Ephonon, m)) * a
    
def impNew(k):
    return k * (numpy.random.randint(0, 2) * 2.0 - 1.0)

def Ek(k, m):
    return hbar ** 2 * k ** 2 * q / (2 * m)

def K(E, m):
    return numpy.sqrt(2 * m * E / (hbar ** 2 * q)) * (numpy.random.randint(0, 2) * 2.0 - 1.0)

def D(E, m):
    return 2 * numpy.sqrt(E * q) / (numpy.pi * (hbar * q) ** 2)

def d_scatter(k, m, n=1e19):
    # n will need to be changed once we have messing to calculate the local carrier densities
    k = abs(k)
    kd = numpy.sqrt(q * n / (epsilon0 * kb * Temp))
    return 2 * ni * q ** 2 * k / (numpy.pi * hbar ** 3 * epsilon0 * kd ** 2 * (kd ** 2 + 4 * k ** 2))

sigmaEp = 18.2e-3
def EphotonDist():
    Epb = -1
    while Epb < 0:
        Epb = numpy.random.normal(Ephoton, sigmaEp, 1) - Eg
    return Epb

def gen(car):
    car.x = abs(numpy.random.normal(0, sigmax, 1))
    car.t0 = numpy.random.normal(0, sigmat, 1)
    car.E = EphotonDist()
    car.kx = K(car.E, me)
    return

def addCarriers(T):
    while carrierReserve[-1].t0[0] < T:
        e = carrierReserve.pop()
        h = copy.deepcopy(e)
        h.charge = 1 * q
        h.m = mh
        h.E = Ek(h.kx, mh)
        carrierList.append(e)
        carrierList.append(h)
        if len(carrierReserve) < 1:
            return 
    return
    
def proploop(carr, Dt=Deltat):
    if carr.tf < Dt:
        carr.propergate(carr.tf)
        carr.drift(carr.tf)
        carr.scatter()
        Dt -= carr.tf
        carr.freeflight()
        if carr.tf < Dt:
            proploop(carr, Dt)
    carr.tf -= Deltat
    carr.propergate(Deltat)
    carr.drift(Deltat)
    return

def makeMesh(carrierList):
    mesh[:] = meshEmpty[:]
    for car in carrierList:
        x = car.x * 1e7 + 1000
        if abs(x) >= len(mesh):
            continue
        mesh[int(x)] += car.charge
    return

def makePotential(mesh, pmesh):
    pmesh[:] = meshEmpty[:]
    for i, Q in enumerate(mesh):
        r = abs(numpy.arange(0, len(pmesh)) - i)
        P = mesh / r
        P[i] = 0
        pmesh[i] = sum(P)
    pmesh *= ke
    return
    
def solveEfield(carrierList):
    makeMesh(carrierList)
    makePotential(mesh, pmesh)
    E = -numpy.gradient(pmesh, dx)
    for car in carrierList:
        x = car.x * 1e7 + 1000
        if abs(x) >= len(E):
            continue
        car.Ex = E[int(x)]
    return

def test(carrierList, T):
    ktot = 0
    kTote = 0
    kToth = 0
    for car in carrierList:
        if car.charge < 0:
            ktot += abs(car.kx)
            kTote += car.kx
        else:
            kToth += car.kx
    kToteL.append(kTote)
    kTothL.append(kToth)
    
    if len(carrierList) > 0:
        kxave.append(2 * ktot / len(carrierList))
        Tk.append(T)
    return

def getX(carrierList):
    xh = []
    xe = []
    for car in carrierList:
        if car.charge > 0:
            xh.append(float(car.x))
        else:
            xe.append(float(car.x))
    Xh.append(xe)
    Xe.append(xh)


map(gen, carrierReserve)  # Makes carrier reserve
carrierReserve = sorted(carrierReserve, key=lambda car:-car.t0)  # Put in generation time order

T = min(carrierReserve[-1].t0)

# 1/0
kxave = []
Tk = []

kToteL = []
kTothL = []


KlistTotal = []
KlistHalf = []
Xh = []
Xe = []
meshArray = []

while T < Tmax:
    if len(carrierReserve) > 0:
        addCarriers(T)
    solveEfield(carrierList)
    carrierList = [carr for carr in carrierList if not carr.destroy]  # Remove carriers that have recombined
#    test(carrierList, T)
    getX(carrierList)
    map(proploop, carrierList)
    T += Deltat
    meshArray.append(mesh)
    print Tmax - T
    
# with open('xe1e61umdeftest', mode='a+b') as file:
#     pickle.dump(Xe, file)   
    
# with open('xh1e61umdeftest', mode='a+b') as file:
#     pickle.dump(Xh, file)   

with open('density', mode='a+b') as file:
    pickle.dump(meshArray, file)   


os.system("ssmtp mark.e.barnes@gmail.com < mail.txt")
print 'Done'
