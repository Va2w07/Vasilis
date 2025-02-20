import pickle
import numpy
import pylab

with open('teflon_elena.pic') as f1:
    s1n = pickle.load(f1)
    s1pe = pickle.load(f1)
    s1ae = pickle.load(f1)
    s1aei = pickle.load(f1)
    s1f = pickle.load(f1)
    s1rat = pickle.load(f1)

with open('teflonreal.pic') as f1:
    s2n = pickle.load(f1)
    s2pe = pickle.load(f1)
    s2ae = pickle.load(f1)
    s2aei = pickle.load(f1)
    s2f = pickle.load(f1)
    s2rat = pickle.load(f1)

# with open('HDPE3.pic') as f1:
#     s3n = pickle.load(f1)
#     s3pe = pickle.load(f1)
#     s3ae = pickle.load(f1)
#     s3aei = pickle.load(f1)
#     s3f = pickle.load(f1)
#     s3rat = pickle.load(f1)

# with open('HDPE12.pic') as f1:
#     s12n = pickle.load(f1)
#     s12pe = pickle.load(f1)
#     s12ae = pickle.load(f1)
#     s12aei = pickle.load(f1)
#     s12f = pickle.load(f1)
#     s12rat = pickle.load(f1)

# with open('HDPE22.pic') as f1:
#     s22n = pickle.load(f1)
#     s22pe = pickle.load(f1)
#     s22ae = pickle.load(f1)
#     s22aei = pickle.load(f1)
#     s22f = pickle.load(f1)
#     s22rat = pickle.load(f1)

# with open('HDPE32.pic') as f1:
#     s32n = pickle.load(f1)
#     s32pe = pickle.load(f1)
#     s32ae = pickle.load(f1)
#     s32aei = pickle.load(f1)
#     s32f = pickle.load(f1)
#     s32rat = pickle.load(f1)

s1f = numpy.array(s1f[0:len(s1f)/2])
s1n = numpy.array(s1n[0:len(s1n)/2])
s1aei = s1aei[0:len(s1aei)/2]
s1ae = s1ae[0:len(s1ae)/2]
s1rat = s1rat[0:len(s1rat)/2]

s2f = numpy.array(s2f[0:len(s2f)/2])
s2n = numpy.array(s2n[0:len(s2n)/2])
s2aei = s2aei[0:len(s2aei)/2]
s2ae = s2ae[0:len(s2ae)/2]
s2rat = s2rat[0:len(s2rat)/2]

# s3f = numpy.array(s3f[0:len(s3f)/2])
# s3n = numpy.array(s3n[0:len(s3n)/2])
# s3aei = s3aei[0:len(s3aei)/2]
# s3ae = s3ae[0:len(s3ae)/2]
# s3rat = s3rat[0:len(s3rat)/2]

# s12f = numpy.array(s12f[0:len(s12f)/2])
# s12n = numpy.array(s12n[0:len(s12n)/2])
# s12aei = s12aei[0:len(s12aei)/2]
# s12ae = s12ae[0:len(s12ae)/2]
# s12rat = s12rat[0:len(s12rat)/2]

# s22f = numpy.array(s22f[0:len(s22f)/2])
# s22n = numpy.array(s22n[0:len(s22n)/2])
# s22aei = s22aei[0:len(s22aei)/2]
# s22ae = s22ae[0:len(s22ae)/2]
# s22rat = s22rat[0:len(s22rat)/2]

# s32f = numpy.array(s32f[0:len(s32f)/2])
# s32n = numpy.array(s32n[0:len(s32n)/2])
# s32aei = s32aei[0:len(s32aei)/2]
# s32ae = s32ae[0:len(s32ae)/2]
# s32rat = s32rat[0:len(s32rat)/2]

L = 21e-3

dx1 = L * (1 - (1 * numpy.cos(s1pe) ) / numpy.sqrt(s1n**2 - 1**2*numpy.sin(s1pe)**2) )
dx2 = L * (1 - (1 * numpy.cos(s2pe) ) / numpy.sqrt(s2n**2 - 1**2*numpy.sin(s2pe)**2) )
# dx3 = L * (1 - (1 * numpy.cos(s3pe) ) / numpy.sqrt(s3n**2 - 1**2*numpy.sin(s3pe)**2) )

# dx12 = L * (1 - (1 * numpy.cos(s12pe) ) / numpy.sqrt(s12n**2 - 1**2*numpy.sin(s12pe)**2) )
# dx22 = L * (1 - (1 * numpy.cos(s22pe) ) / numpy.sqrt(s22n**2 - 1**2*numpy.sin(s22pe)**2) )
# dx32 = L * (1 - (1 * numpy.cos(s32pe) ) / numpy.sqrt(s32n**2 - 1**2*numpy.sin(s32pe)**2) )

A1 = numpy.sqrt(1+(dx1*numpy.pi*s1f*1e12*s1ae/(4*3e8))**2)
A2 = numpy.sqrt(1+(dx2*numpy.pi*s2f*1e12*s2ae/(4*3e8))**2)
# A3 = numpy.sqrt(1+(dx3*numpy.pi*s3f*1e12*s3ae/(4*3e8))**2)

# A12 = numpy.sqrt(1+(dx12*numpy.pi*s12f*1e12*s12ae/(4*3e8))**2)
# A22 = numpy.sqrt(1+(dx22*numpy.pi*s22f*1e12*s22ae/(4*3e8))**2)
# A32 = numpy.sqrt(1+(dx32*numpy.pi*s32f*1e12*s32ae/(4*3e8))**2)

theta21 = numpy.arcsin(numpy.sin(s1ae)*1/s1n)
theta22 = numpy.arcsin(numpy.sin(s2ae)*1/s2n)
# theta23= numpy.arcsin(numpy.sin(s3ae)*1/s3n)

# theta212 = numpy.arcsin(numpy.sin(s12ae)*1/s12n)
# theta222 = numpy.arcsin(numpy.sin(s22ae)*1/s22n)
# theta232= numpy.arcsin(numpy.sin(s32ae)*1/s32n)

T11 = 2*numpy.cos(s1ae)/(numpy.cos(theta21)+s1n*numpy.cos(s1ae))
T21 = 2*s1n*numpy.cos(theta21)/(s1n*numpy.cos(s1ae)+numpy.cos(theta21))
T12 = 2*numpy.cos(s2ae)/(numpy.cos(theta22)+s2n*numpy.cos(s2ae))
T22 = 2*s2n*numpy.cos(theta22)/(s2n*numpy.cos(s2ae)+numpy.cos(theta22))
# T13 = 2*numpy.cos(s3ae)/(numpy.cos(theta23)+s3n*numpy.cos(s3ae))
# T23 = 2*s3n*numpy.cos(theta23)/(s3n*numpy.cos(s3ae)+numpy.cos(theta23))

# T112 = 2*numpy.cos(s12ae)/(numpy.cos(theta212)+s12n*numpy.cos(s12ae))
# T212 = 2*s12n*numpy.cos(theta212)/(s12n*numpy.cos(s12ae)+numpy.cos(theta212))
# T122 = 2*numpy.cos(s22ae)/(numpy.cos(theta222)+s22n*numpy.cos(s22ae))
# T222 = 2*s22n*numpy.cos(theta222)/(s22n*numpy.cos(s22ae)+numpy.cos(theta222))
# T132 = 2*numpy.cos(s32ae)/(numpy.cos(theta232)+s32n*numpy.cos(s32ae))
# T232 = 2*s32n*numpy.cos(theta232)/(s32n*numpy.cos(s32ae)+numpy.cos(theta232))

TT1 = T11*T21
TT2 = T12*T22
# TT3 = T13*T23

# TT12 = T112*T212
# TT22 = T122*T222
# TT32 = T132*T232

w1 = 2*3e8/(s1f*1e12*numpy.pi*s1ae)*numpy.sqrt(1+(dx1*s1f*numpy.pi*s1ae**2/(4*3e8))**2)
w2 = 2*3e8/(s2f*1e12*numpy.pi*s2ae)*numpy.sqrt(1+(dx2*s2f*numpy.pi*s2ae**2/(4*3e8))**2)
# w3 = 2*3e8/(s3f*1e12*numpy.pi*s3ae)*numpy.sqrt(1+(dx*s3f*numpy.pi*s3ae**2/(4*3e8))**2)

# w12 = 2*3e8/(s12f*1e12*numpy.pi*s12ae)*numpy.sqrt(1+(dx12*s12f*numpy.pi*s12ae**2/(4*3e8))**2)
# w22 = 2*3e8/(s22f*1e12*numpy.pi*s22ae)*numpy.sqrt(1+(dx*s22f*numpy.pi*s22ae**2/(4*3e8))**2)
# w32 = 2*3e8/(s32f*1e12*numpy.pi*s32ae)*numpy.sqrt(1+(dx*s32f*numpy.pi*s32ae**2/(4*3e8))**2)

pylab.figure()
pylab.plot(s1f, TT1/A1, label = 'E with T' )
pylab.plot(s1f, 1.0/A1, label = 'E' )
pylab.plot(s1f, numpy.abs(s1rat), label = 'E Amplitude of T ratio')
pylab.ylim([0.0,2])
pylab.xlim([0.1,2.5])
pylab.legend()

pylab.figure()
pylab.plot(s2f, TT2/A2, label = 'S With T' )
pylab.plot(s2f, 1.0/A2, label = 'S' )
pylab.plot(s2f, numpy.abs(s2rat), label = 'S Amplitude of T ratio')
pylab.xlim([0.1,2.5])
pylab.ylim([0.0,2])
pylab.legend()

# pylab.figure()
# pylab.plot(s3f, TT3/A3, label = 'D3 Eq. 10 With T' )
# pylab.plot(s3f, 1.0/A3, label = 'D3 Eq. 10' )
# pylab.plot(s3f, numpy.abs(s3rat), label = 'D3 Amplitude of T ratio')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,2])
# pylab.legend()


# pylab.figure()
# #pylab.plot(s12f, TT12/A12, label = 'D1 Eq. 10 With T 2' )
# pylab.plot(s12f, 1.0/A12, label = 'D1 Eq. 10 2' )
# pylab.plot(s12f, numpy.abs(s12rat), label = 'D1 Amplitude of T ratio 2')
# pylab.ylim([0.0,2])
# pylab.xlim([0.5,2.5])
# pylab.legend()


# pylab.figure()
# #pylab.plot(s22f, TT12/A22, label = 'D2 Eq. 10 With T 2' )
# pylab.plot(s22f, 1.0/A22, label = 'D2 Eq. 10 2' )
# pylab.plot(s22f, numpy.abs(s22rat), label = 'D2 Amplitude of T ratio 2')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,2])
# pylab.legend()


# pylab.figure()
# #pylab.plot(s32f, TT12/A32, label = 'D3 Eq. 10 With T 2' )
# pylab.plot(s32f, 1.0/A32, label = 'D3 Eq. 10 2' )
# pylab.plot(s32f, numpy.abs(s32rat), label = 'D3 Amplitude of T ratio 2')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,2])
# pylab.legend()

# pylab.figure()
# pylab.plot(s1f, A * 2*3e8/(s1f*numpy.pi*s1ae*1e12)  )

# pylab.figure()
# pylab.plot(s1f, w)

pylab.figure()
pylab.plot(s1f, numpy.abs(s1rat), s1f, numpy.angle(s1rat) )

pylab.figure()
pylab.plot(s1f,s1n.real, s1f, s1n.imag, label = 'Elena')
pylab.plot(s2f,s2n.real, s2f, s2n.imag, label = 'Sam')
# pylab.plot(s3f,s3n.real, s3f, s3n.imag)
pylab.xlim([0.1,3.0])

pylab.figure()
pylab.plot(s1f, s1ae.real, label = 'Elena Real')
pylab.plot(s2f, s2ae.real, label = 'Sam Real')
# pylab.plot(s3f, s3ae.real, label = 'Data 3, Focus 1')
pylab.xlim([0.1,2.5])
pylab.ylim([0.0,0.4])

pylab.legend()

# pylab.figure()
# pylab.plot(s12f, s12ae.real, label = 'Data 1, Focus 2')
# pylab.plot(s22f, s22ae.real, label = 'Data 2, Focus 2')
# pylab.plot(s32f, s32ae.real, label = 'Data 3, Focus 2')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,0.2])
# pylab.legend()

# pylab.figure()
# pylab.plot(s1f, s1ae.real, label = 'Data 1, Focus 1')
# pylab.plot(s2f, s2ae.real, label = 'Data 2, Focus 1')
# pylab.plot(s3f, s3ae.real, label = 'Data 3, Focus 1')

# pylab.plot(s12f, s12ae.real, label = 'Data 1, Focus 2')
# pylab.plot(s22f, s22ae.real, label = 'Data 2, Focus 2')
# pylab.plot(s32f, s32ae.real, label = 'Data 3, Focus 2')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,0.2])
# pylab.legend()


# pylab.figure()
# pylab.plot(s1f, s1ae.real, label = 'Data 1, Focus 1')
# pylab.plot(s12f, s12ae.real, label = 'Data 1, Focus 2')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,0.2])
# pylab.legend()

# pylab.figure()
# pylab.plot(s2f, s2ae.real, label = 'Data 2, Focus 1')
# pylab.plot(s22f, s22ae.real, label = 'Data 2, Focus 2')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,0.2])
# pylab.legend()

# pylab.figure()
# pylab.plot(s3f, s3ae.real, label = 'Data 3, Focus 1')
# pylab.plot(s32f, s32ae.real, label = 'Data 3, Focus 2')
# pylab.xlim([0.5,2.5])
# pylab.ylim([0.0,0.2])
# pylab.legend()



# pylab.figure()
# pylab.plot(s1f, s1aei.real)


pylab.figure()
pylab.plot(s1f, s1ae.imag, label = 'new')
pylab.plot(s1f, s1aei.imag, label = 'old')
pylab.plot(s1f, s1pe.imag, label = 'Aaron')
pylab.xlim([0.5,2.5])

pylab.legend()

pylab.figure()
pylab.plot(s2f, s2ae.imag, label = 'new')
pylab.plot(s2f, s2aei.imag, label = 'old')
pylab.plot(s2f, s2pe.imag, label = 'Aaron')
pylab.xlim([0.5,2.5])

pylab.legend()

pylab.show()
