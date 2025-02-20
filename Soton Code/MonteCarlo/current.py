import numpy
import pylab
import pickle 

dt = 1e-14
xe = []
xh = []

L = 30e-6
dx = 10e-9

earray = []
harray = []


with open( "xe1e61um10fs", "rb" ) as file:
    while 1:
        try:
            xe.append(pickle.load(file))
        except EOFError:
            break

with open( "xh1e61um10fs", "rb" ) as file:
    while 1:
        try:
            xh.append(pickle.load(file))
        except EOFError:
            break
i = 0

for elist in xe[-1]:
    ebin = numpy.zeros(3000)
    for e in elist:
        if e != []:
            e = e * 1e8
            e = e + 1000
            ebin[int(e)] += 1
    print i
    i += 1
    earray.append(ebin)

i = 0
for hlist in xh[-1]:
    hbin = numpy.zeros(3000)
    for h in hlist:
        if h != []:
            h = h*1e8
            h = h+1000
            hbin[int(h)] += 1
    print i
    i += 1
    harray.append(hbin)

earray = numpy.array(earray)
harray = numpy.array(harray)

density = earray-harray

all = []
half = []
for den in density:
    all.append(sum(den))
    half.append(sum(den[1000:len(den)-1]))

X = numpy.arange(0,dx*len(density[7]), dx)*1e6-10
pylab.plot(X, density[7])
#pylab.plot(numpy.gradient(numpy.gradient(numpy.array(all))))
#pylab.plot(numpy.gradient(numpy.gradient(numpy.array(half))))
pylab.savefig('current_density', format='pdf')