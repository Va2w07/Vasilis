import pickle
import numpy
import pylab

thetarange = numpy.arange(0, 1, 0.01)
n = []
pe = []
ae = []
f = []
diffang = []
diffarray = []
aang = []
pang = []
fpos = 439

for theta in thetarange:
    with open('ang2'+str(theta)+'.pic') as f1:
        n.append(pickle.load(f1))
        pe.append(pickle.load(f1))
        aetemp = pickle.load(f1)
        f.append(pickle.load(f1))
        aetemp = aetemp[0:len(pe[-1])] 
        ae.append(aetemp)
        diffang.append(numpy.array(pe[-1]) - numpy.array(ae[-1]))
        aang.append(ae[-1][fpos])
        pang.append(pe[-1][fpos]) 
        diffarray.append(diffang[-1][fpos])

pylab.plot(f[20][0:len(pe[-1])], ae[20], label = 'Analytical')
pylab.plot(f[20][0:len(pe[-1])], pe[20], label = 'Aaron')
pylab.plot(f[20][0:len(pe[-1])], diffang[20], label = 'Diff')
pylab.legend()
pylab.show()
pylab.plot(thetarange, aang, label = 'Analytical')
pylab.plot(thetarange, pang, label = 'Aaron')
pylab.legend()
pylab.show()

pylab.plot(thetarange, abs(numpy.array(diffarray)))
pylab.show()

