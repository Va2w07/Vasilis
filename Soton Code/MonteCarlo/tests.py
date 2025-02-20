import numpy
import pylab
import pickle 

dt = 1e-14
xe = []
xh = []

with open( "xe1e61umphonon", "rb" ) as file:
    while 1:
        try:
            xe.append(pickle.load(file))
        except EOFError:
            break

with open( "xh1e61umphonon", "rb" ) as file:
    while 1:
        try:
            xh.append(pickle.load(file))
        except EOFError:
            break

ce = []
cea = []
qe = []
qea = []
for elist in xe[-1]:
    if len(elist) > 0:
        print len(elist)
        qe.append(sum(elist))
        ce.append(qe[-1]/len(elist))
    ae = [e for e in elist if e > 0]
    if len(ae) > 0:
        qea.append(sum(ae))
        cea.append(qea[-1]/len(ae))
        
ch = []
cha = [] 
qh = []
qha = []

for hlist in xh[-1]:
    if len(hlist) > 0:
        qh.append(sum(hlist))
        ch.append(qh[-1]/len(hlist))

    ah = [h for h in hlist if h > 0]
    if len(ah) > 0:
        qha.append(sum(ah))
        cha.append(qha[-1]/len(ah))

cha = numpy.array(cha)
cea = numpy.array(cea)

ch = numpy.array(ch)
ce = numpy.array(ce)

qe = numpy.array(qe)
qea = numpy.array(qea)
qh = numpy.array(qh)
qha = numpy.array(qha)

print len(cha), len(cea)
#1/0
da = cha-cea
d = ch-ce

qa = qha+qea
q = qh+qe

P = q*d
Pa = qa*da

T = numpy.arange(0,dt*(len(P)),dt)
print len(T), len(P)
pylab.plot(T,numpy.gradient(numpy.gradient(P)))
pylab.plot(T,numpy.gradient(numpy.gradient(Pa)))
#pylab.savefig('THz_MC', format='pdf')
pylab.show()