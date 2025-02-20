import extraction
import csv
import pylab
import numpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt 
def quartz(ref,samp,thickness):
	# create two arrays tArray 0 and 1 with the data of the ref file
	tArray = [list(),list(),list(),list()]
	with open(ref) as csvfile:
		csvreader =  csv.reader(csvfile, delimiter = ',')
		for row in csvreader:
			tArray[0].append(float(row[0]))
			tArray[1].append(float(row[1]))
	
	# create two arrays tArray 2 and 3 with the data of the sample file		
	with open(samp) as csvfile:
		csvreader =  csv.reader(csvfile, delimiter = ',')
		for row in csvreader:
			tArray[2].append(float(row[0]))
			tArray[3].append(float(row[1]))
	print tArray[0]
	print tArray[1]
	print tArray[2]
	print tArray[3]

	exp = int(numpy.log2(len(tArray[0])))+5
	L = 2**exp
	dt=((tArray[0][-1]-tArray[0][0])*1e-12/(len(tArray[0])))#average dt
	print dt
	xmaxr=numpy.argmax(tArray[1])
	xmaxs=numpy.argmax(tArray[3])
	offset = numpy.abs(int((tArray[0][0] - tArray[2][0])*1e-12/dt))
	tch = ((numpy.abs(((tArray[0][0] - tArray[2][0])*1e-12/dt)))-(numpy.abs(int((tArray[0][0] - tArray[2][0])*1e-12/dt))))
	#dt=(((tArray[2][1]-tArray[2][0])+ (tArray[0][1]-tArray[0][0]))/2)*1e-12
	#print offset
	#print tch 
	
	#dt = 19.0e-15
	#dt2=(tArray[2][-1]-tArray[2][0])*1e-12/len(tArray[2])
	#dt = (tArray[0][1]-tArray[0][0])*1e-12
	#print dt
    
	ddt=tch*dt
	print ddt
	t1 = numpy.arange(0, L, 1) * dt
	t2 = numpy.arange(offset, L+offset, 1) * dt
	wr = window(xmaxr,tArray[1])
	ws = window(xmaxs,tArray[3])
	ref = (numpy.array(tArray[1]) - tArray[1][0])*wr
	sample = (numpy.array(tArray[3])  - tArray[3][0])*ws

	fref = interp1d((numpy.array(tArray[0])-tArray[0][0])*1e-12, ref)
	fsam = interp1d((numpy.array(tArray[2]) - tArray[2][0])*1e-12, sample)
	# fref = (numpy.array(tArray[0])-tArray[0][0])*1e-12
	# fsam = (numpy.array(tArray[2]) - tArray[2][0])*1e-12
	pylab.plot(tArray[0], ref, tArray[2], sample)

	ref = fref(numpy.arange(0,len(ref))*dt)
	sample = fsam(numpy.arange(0,len(sample))*dt)
	(xRVal2, yRVal2, xSVal2, ySVal2, nReal, nImag, f, start, end) = extraction.getRefrac(t1, ref, t2, sample, thickness, 1,ddt,dt, L)
	# pylab.plot(f*1e-12, nReal[0:4096], f*1e-12, nImag[0:4096])
	# pylab.plot(tArray[0], tArray[1])
	pylab.show()
	return (nReal, nImag, f)

def window(xmax,x):
	width = 150.0
	sigma = width/numpy.sqrt(2*numpy.log(2))
	w=numpy.ones(len(x))
	for x in range (0,len(w)):
		if (x > (xmax+20)):
			w[x]=numpy.exp((-(x-xmax)**2)/(2*sigma**2))
	return w





n, ni,fx = quartz('2mmr_1.csv','2mms.txt', 2.000)

plt.plot(fx*1e-12,n, linewidth=2.0, label='real ref index', color='red', linestyle='-')
plt.show()
plt.plot(fx*1e-12,ni, linewidth=2.0, label='imag ref index', color='blue', linestyle='--')



plt.show()


# nh, nih,fh = quartz('slabhept1.txt','hept1.txt', 2.292)

# plt.plot(fh,nh, linewidth=2.0, label='real ref index', color='red', linestyle='-')
# plt.show()
# plt.plot(fh,nih, linewidth=2.0, label='imag ref index', color='blue', linestyle='--')



# plt.show()











