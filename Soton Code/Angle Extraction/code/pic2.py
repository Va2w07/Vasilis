import pickle
import numpy
import pylab

with open('cretes1.pic') as f1:
    s1n = pickle.load(f1)
    s1pe = pickle.load(f1)
    s1ae = pickle.load(f1)
    s1f = pickle.load(f1)
s1f = s1f[0:len(s1f)/2]
s1ae = s1ae[0:len(s1ae)/2]

with open('cretes2.pic') as f1:
    s2n = pickle.load(f1)
    s2pe = pickle.load(f1)
    s2ae = pickle.load(f1)
    s2f = pickle.load(f1)
s2f = s2f[0:len(s2f)/2]
s2ae = s2ae[0:len(s2ae)/2]

with open('cretes3.pic') as f1:
    s3n = pickle.load(f1)
    s3pe = pickle.load(f1)
    s3ae = pickle.load(f1)
    s3f = pickle.load(f1)
s3f = s3f[0:len(s3f)/2]
s3ae = s3ae[0:len(s3ae)/2]

with open('cretep1.pic') as f1:
    p1n = pickle.load(f1)
    p1pe = pickle.load(f1)
    p1ae = pickle.load(f1)
    p1f = pickle.load(f1)
p1f = p1f[0:len(p1f)/2]
p1ae = p1ae[0:len(p1ae)/2]

with open('cretep2.pic') as f1:
    p2n = pickle.load(f1)
    p2pe = pickle.load(f1)
    p2ae = pickle.load(f1)
    p2f = pickle.load(f1)
p2f = p2f[0:len(p2f)/2]
p2ae = p2ae[0:len(p2ae)/2]

with open('cretep3.pic') as f1:
    p3n = pickle.load(f1)
    p3pe = pickle.load(f1)
    p3ae = pickle.load(f1)
    p3f = pickle.load(f1)
p3f = p3f[0:len(p3f)/2]
p3ae = p3ae[0:len(p3ae)/2]

with open('axel.pic') as f1:
    an = pickle.load(f1)
    ape = pickle.load(f1)
    aae = pickle.load(f1)
    af = pickle.load(f1)
af = af[0:len(af)/2]
aae = aae[0:len(aae)/2]

astart = 50
aend= 200

cstart = 50
csend = 200

#s1f = s1f[cstart:csend]
#s2f = s2f[cstart:csend]
#s3f = s3f[cstart:csend]

#p1f = p1f[cstart:csend]
#p2f = p2f[cstart:csend]
#p3f = p3f[cstart:csend]

#af = af[cstart:csend]

#s1pe = s1pe[cstart:csend]
#s2pe = s2pe[cstart:csend]
#s3pe = s3pe[cstart:csend]

#p1pe = p1pe[cstart:csend]
#p2pe = p2pe[cstart:csend]
#p3pe = p3pe[cstart:csend]

#ape = ape[cstart:csend]

#s1ae = s1ae[cstart:csend]
#s2ae = s2ae[cstart:csend]
#s3ae = s3ae[cstart:csend]

#p1ae = p1ae[cstart:csend]
#p2ae = p2ae[cstart:csend]
#p3ae = p3ae[cstart:csend]

#aae = aae[cstart:csend]

#s1n = s1n[cstart:csend]
#s2n = s2n[cstart:csend]
#s3n = s3n[cstart:csend]

#p1n = p1n[cstart:csend]
#p2n = p2n[cstart:csend]
#p3n = p3n[cstart:csend]

#an = an[cstart:csend]


pylab.plot(s1f, s1pe, label = 's1')
pylab.plot(s2f, s2pe, label = 's2')
pylab.plot(s3f, s3pe, label = 's3')

pylab.plot(p1f, p1pe, label = 'p1')
pylab.plot(p2f, p2pe, label = 'p2')
pylab.plot(p3f, p3pe, label = 'p3')
pylab.legend()
pylab.show()

pylab.plot(s1f, s1pe, label = 's1 aaron')
pylab.plot(s1f, s1ae, label = 's1 analytical')
pylab.legend()
pylab.show()
print len(af), len(ape)
pylab.plot(af, ape, label = 'axel arron')
pylab.plot(af, aae, label = 'axel analytical')
pylab.legend()
pylab.show()

pylab.plot(s1f, s1pe - s1ae, label = 's1')
pylab.plot(s2f, s2pe - s2ae, label = 's2')
pylab.plot(s3f, s3pe - s3ae, label = 's3')

pylab.plot(p1f, p1pe - p1ae, label = 'p1')
pylab.plot(p2f, p2pe - p2ae, label = 'p2')
pylab.plot(p3f, p3pe - p3ae, label = 'p3')
pylab.legend()
pylab.show()

pylab.plot(af, ape - aae, label = 'axel diff')
pylab.legend()
pylab.show()

pylab.plot(s1pe, s1pe - s1ae,'o')
pylab.show()
#pylab.plot(s2pe, s2pe - s2ae,'o')
#pylab.plot(s3pe, s3pe - s3ae,'o')

#pylab.plot(p1pe, p1pe - p1ae,'o')
#pylab.plot(p2pe, p2pe - p2ae,'o')
#pylab.plot(p3pe, p3pe - p3ae,'o')

pylab.plot(ape, ape - aae,'o')
pylab.show()
