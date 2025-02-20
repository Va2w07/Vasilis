import numpy
import pylab
import param
from numpy.fft import fft

def slab1(n):
    return param.logT1(1.0, n, param.kD[0]) - param.T0

def slab1dash(n):
    return param.logT1dash(1.0, n, param.kD[0])


def bufferData(xR, yR, xS, yS):
    
    dtR = (xR[-1] - xR[0])/(len(xR)) 
    dtC = (xS[-1] - xS[0])/(len(xS)) 
      
    minT = min(xR[0], xS[0])
    maxT = max(xR[-1], xS[-1]) 
        
    paddingR = (int((xR[0] - minT)/dtR),int((maxT - xR[-1])/dtR))
    paddingC = (int((xS[0] - minT)/dtC),int((maxT - xS[-1])/dtC))

    yR = numpy.append(yR,0)    
    yS = numpy.append(yS,0)    
    
    yS = numpy.insert(yS,0,0)    
    yR = numpy.insert(yR,0,0)    
    
    yR = numpy.lib.pad(yR,paddingR,'edge')
    yS = numpy.lib.pad(yS,paddingC,'edge')

    M = int(numpy.log2(len(yS)))+3
    extra = 2**M - len(yS)

    yS = numpy.lib.pad(yS,(0,extra),'edge')
    yR = numpy.lib.pad(yR,(0,extra),'edge')

    return yR,yS

def getGoodSignal(kstart, YS):

    weight=numpy.zeros(len(YS), float)
    YSMag = abs(YS)
    peak = YSMag[kstart]
    minMag = peak*0.1
    
    start = 0
    end = 0
    
    for i in range(kstart-1, 0,-1):
        if YSMag[i] < minMag:
            start = i
            break

    for i in range(kstart, len(YS)/2):
        if YSMag[i] < minMag:
            end = i
            break
    
    sigma1 = (kstart - start)/(numpy.sqrt(2*numpy.log(10)))
    sigma2 = (end - kstart)/(numpy.sqrt(2*numpy.log(10)))
             
    for k in range(0,len(weight)/2):
        if k < kstart:
            c = sigma1
        else:
            c = sigma2
        weight[k]=numpy.exp(-(k-kstart)**2/(2*float(c)**2))
        weight[len(YS)-k-1]=weight[k]

    return weight, start, end


def getRefrac(xR, yR, xS, yS, thickness, tUnit):

    thickness=numpy.array([thickness*1e-3])
    
    xR = numpy.array(xR)*1e-12
    xS = numpy.array(xS)*1e-12
    yR = numpy.array(yR)
    yS = numpy.array(yS)
        
    maxy = max(max(yR),max(yS))
    
    yR = yR/maxy
    yS = yS/maxy

    dt = (xR[-1] - xR[0])/len(xR)   

    lenOrginal = len(xR) 
   
    (yR, yS) = bufferData(xR,yR,xS,yS)
    
    N = len(yR)
    M = N/2.0
       
    df=1.0/(N*dt)
    dk = 2*numpy.pi*df/2.9979e8
    
    time = numpy.arange(0,N)*dt
    f = numpy.arange(0,M)*df
    
    xR = time
    xS = time
   
    YR = fft(yR)/N
    YS = fft(yS)/N
    W = YS/YR
    
    logW=numpy.log(W[0:M])
   
    kstart = numpy.argmax(numpy.abs(YS[1:M]))
    
    (weight, start, end) = getGoodSignal(kstart, YS)
       
    logWfixed=param.pahseUnqwarp(logW, kstart, weight, M)
    
    (k_lo, k_hi, n,tfdiff, tfexp) = param.convert_to_ri(kstart, thickness, slab1, slab1dash, logWfixed, dk)
    
    return xR[0:lenOrginal]*1e12, yR[0:lenOrginal], xS[0:lenOrginal]*1e12, yS[0:lenOrginal], n.real, n.imag, f*1e-12, start, end

    
