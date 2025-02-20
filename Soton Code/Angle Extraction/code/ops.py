import numpy
import FileDialog
import pylab
import csv
import models

# remove fft defined in physics.py

from numpy.fft import fft, ifft
#from cmath import *      # note this prevents Ufuncs existing
#
global T0, kD, freqk
T0=0.0; kD=[0.0]

#------------------------------------------------------------------------------
# Read the time base of the data file and determine if the sample file needs
# to be shifted.
#--------------------------------------------------------------------------


def findtimeshift(reference_file, sample_file):
    'import the time base of sample and reference file'
    tx=readdatatimebase(reference_file, 1.0)
    ty=readdatatimebase(sample_file, 1.0)
    'determine the time shift'
    i=0
    while tx[i]<ty[0]:
        i+=1
        #print tx[i],ty[0]
    'determine which element is closest'
    if abs(ty[0]-tx[i-1])<abs(ty[0]-tx[i]):
        i-=1
    return i


def readdatafile(filename, scale, paddingfactor):
    "reads datafile and appends same length of zeros"
    z=[]
    fileobj=open(filename, 'r')
    for line in fileobj.readlines():
        z.append(scale*float(line.split(',')[1]))
    z = numpy.array(z)
    z = z-z[0]
    z = z.tolist()
    for i in range(0,paddingfactor*len(z)):
        z.append(0.0)
    return numpy.array(z)

def readdatafile2(filename, scale):
    "reads datafile with no padding"
    z=[]
    fileobj=open(filename, 'r')
    for line in fileobj.readlines():
        z.append(scale*float(line.split(',')[1]))
    return numpy.array(z)

def readdatatimebase(filename, scale):
    "reads timebase of data file"
    z=[]
    fileobj=open(filename, 'r')
    for line in fileobj.readlines():
        z.append(scale*float(line.split(',')[0]))
    for i in range(0,len(z)):
        z.append(0.0)
    return numpy.array(z)

def readdatafilewidth(filename, scale):
    "reads datafile and appends same length of zeros"
    z=[]
    fileobj=open(filename, 'r')
    for line in fileobj.readlines():
        z.append(line.split(',')[1])
    return numpy.array(z)

def savedatafile(y, filename, scale):
    M=len(y)/2
    fileobj=open(filename, 'w')
    for i in range(0,M):
        fileobj.write(str(0) + "    " + str(y[i]*scale) + "\n")
    fileobj.close()
    
def savedatafilefreq(y, filename, paddingfactor, deltaf, multifactor):
    M=len(y)/paddingfactor
    fileobj=open(filename, 'w')
    for i in range(0,M):
        fileobj.write(str(i*deltaf*1e-12*paddingfactor)+","+str(multifactor*y[i*paddingfactor])+"\n")
    fileobj.close() 
        
def savedatafiletime(y, filename, deltat):
    M=len(y)
    fileobj=open(filename, 'w')
    for i in range(0,M):
        fileobj.write(str(i*deltat*1e12)+","+str(y[i])+"\n")
    fileobj.close()  

def truncate(z, length, fall):
    "Removes end of data with smooth cut-off"
    N=len(z)
# find peak
    place=numpy.argmax(abs(z))
    #place = 612
    "Geoffs old code"
    #peak=-1e10; place=0
    #for t in range(0,N):
    #    if abs(z[t])>peak: peak=abs(z[t]); place=t

# set weight function
    "crete_data4"
    weight=time_window(length+place, fall, N)*time_window2(place-length, fall, N)
    "Use this for a prepulse, i.e Cambridge data si, 3si, topas"
    #weight=time_window(length+place, fall, N)*time_window2(place-length, fall, N)
    #weight=time_window2(place-length, fall, N)
    "Use this for normal truncation"
    #weight=time_window(length+place, fall, N)
    #file = open('weight3.txt','w')
    #for i in range(0, len(z)):
    #    file.writelines(str(i*deltat*1e12)+","+str(weight[i])+"\n")
    #file.close()
    #pylab.figure()
    #pylab.plot(tscale,weight)
    #pylab.plot(tscale,z/numpy.max(abs(z)))
    #pylab.plot(tscale,weight*z/numpy.max(abs(z)))
    "sharp cutoff"
    ##weight=numpy.zeros(N)
    ##for i in numpy.arange(0,place+length):
    ##    weight[i]=1.0       
    "write truncation function to file" 
    #file = open('weightsmoothsam'+str(600)+'.txt','w')
    #for i in range(0, len(z)):
    #    file.writelines(str(i*deltat*1e12)+","+str(weight[i])+"\n")
    #file.close()     
    return z*weight


def time_window(where, scale, N):
    "model for removing data with long time delays"
    weight=1.0/(1.0 + numpy.exp((numpy.arange(0,N)-where)/float(scale)))
    return weight

def time_window2(where, scale, N):
    "model for removing data with long time delays"
    weight=1.0/(1.0 + numpy.exp(-(numpy.arange(0,N)-where)/float(scale)))
    return weight

def noise_model_time(z):
    "Totally empirical model to get noise estimate"
    N=len(z)
    Z=fft(abs(z))/N
    K=100
    for i in range(1,K):
        f=numpy.exp(-0.001*i*i)
        Z[i]=Z[i]*f
        Z[N-i]=Z[N-i]*f
    Z[K:N-K]=0.0
    z=N*ifft(Z).real
    return 5e-4 + 0.02*z*z

def good_signal(N, kstart, goodrange):
    "Model of range where signal/noise is good, used for fixing phase"
    weight=numpy.zeros(N, float)
    M=N/2
    for k in range(10,M+1):
        weight[k]=numpy.exp(-((k-kstart)/float(goodrange))**2)
        weight[N-k]=weight[k]
    weight[0]=0.0;
    return weight

def good_signal_new(kstart, YS):
    weight=numpy.zeros(len(YS), float)
    YSMag = abs(YS)
    peak = YSMag[kstart]
    minMag = peak*0.5
    
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
    pylab.figure(101)
    pylab.plot(numpy.log(YS)/max(abs(numpy.log(YS))))
    pylab.plot(weight)

    return weight, start, end


class Phase:
    def __init__(self,x):
        self.x = x
    description = "This shape has not been described yet"
    author = "Nobody has claimed to make this shape yet"
    def fixphase(self, kstart, weight, phasetolerance, phase_offset):
        M=len(self.x)
        phase=numpy.array(self.x.imag)
    # start at freq kstart and step down
        lastph=phase[kstart]
        n=0
        for k in range(kstart, -1, -1):
            if phase[k]-lastph >  phasetolerance*numpy.pi: n=n-1
            if phase[k]-lastph < -phasetolerance*numpy.pi: n=n+1
            lastph=phase[k]
            phase[k]=phase[k] + 2*n*numpy.pi
    # start at freq kstart and step up
        lastph=phase[kstart]
        n=0
        for k in range(kstart, M, 1):
            if phase[k]-lastph >  phasetolerance*numpy.pi: n=n-1
            if phase[k]-lastph < -phasetolerance*numpy.pi: n=n+1
            lastph=phase[k]
            phase[k]=phase[k] + 2*n*numpy.pi
    # extrapolate to zero frequency
    
        krange = numpy.array(range(0,M))
        sumw=sum(weight[0:M])
        sumwk=sum(weight[0:M]*krange)
        sumwksq=sum(weight[0:M]*krange*krange)
        sumwph=sum(weight[0:M]*phase)
        sumwphk=sum(weight[0:M]*krange*phase)
        b=(sumwph*sumwksq-sumwphk*sumwk)/(sumw*sumwksq-sumwk*sumwk)

        phase=phase-b
        
        
        return self.x.real + (1J)*(phase + phase_offset)
    
    def newfixphase(self, W, kstart, weight):
        "unwraps phase using phase differences"
        M=len(self.x)
        phase=numpy.zeros(M, float)
        
#        for k in range(kstart, -1, -1):
#            Z=numpy.conjugate(W[k])*W[k+1]
#            phase[k-1] = phase[k] + numpy.arctan2(Z.imag, Z.real)
        logW = numpy.log(W)
        phase[kstart] = logW[kstart].imag
        
        for k in range(kstart, M-1, 1):
            Z=numpy.conjugate(W[k])*W[k+1]
            phase[k+1] = phase[k] + numpy.arctan2(Z.imag, Z.real)
            # extrapolate to zero frequency

        for k in range(kstart, 0, -1):
            Z=numpy.conjugate(W[k-1])*W[k]
            phase[k-1] = phase[k] - numpy.arctan2(Z.imag, Z.real)

        krange = numpy.array(range(0,M))

        sumw=sum(weight[0:M])
        sumwk=sum(weight[0:M]*krange)        
        sumwksq=sum(weight[0:M]*krange*krange)        
        sumwph=sum(weight[0:M]*phase)
        sumwphk=sum(weight[0:M]*krange*phase)

        b=(sumwph*sumwksq-sumwphk*sumwk)/(sumw*sumwksq-sumwk*sumwk)

        p=round(b/(2*numpy.pi))
        phase=phase-2*p*numpy.pi   
        
        return self.x.real + (1J)*(phase )

    
def convert_to_ri_thinfilm(kstart, thickness, model, deriv, logWfixed, dk, nsubstrate, Ls):
    "Fits model to logT to get refractive index"
    global T0, kD, kDs
    M=len(logWfixed)
    n=numpy.zeros(M,complex)#
    tfdiff=numpy.zeros(M,complex)
    tfexp=numpy.zeros(M,complex)
# start at good freq and step down
    z0=1.0-0.005J
    for k in range(kstart, 0, -1):
        T0=logWfixed[k]
        kD=k*dk*thickness
        kDs=k*dk*Ls
        ns=nsubstrate[k]
        (ok,z0,tfd,tfe)=newtonthinfilm(z0, model, deriv, ns)
        if not ok: break
        n[k]=z0
        tfdiff[k]=tfd
        tfexp[k]=tfe
    k_lo =k
# start at kstart and step up
    z0=1.0-0.008J
    for k in range(kstart, M):
        T0=logWfixed[k]
        kD=k*dk*thickness
        kDs=k*dk*Ls
        ns=nsubstrate[k]
        (ok,z0,tfd,tfe)=newtonthinfilm(z0, model, deriv, ns)
        if not ok: break
        n[k]=z0
        tfdiff[k]=tfd
        tfexp[k]=tfe
    k_hi=k
    return (k_lo, k_hi, n, tfdiff, tfexp)

def convert_to_ri(kstart, thickness, model, deriv, logWfixed, dk, numRefl):
    "Fits model to logT to get refractive index"
    global T0, kD
    M=len(logWfixed)
    n=numpy.zeros(M,complex)#
    tfdiff=numpy.zeros(M,complex)
    tfexp=numpy.zeros(M,complex)
# start at good freq and step down
    z0=2.0-0.005J
    for k in range(kstart, 0, -1):
        T0=logWfixed[k]
        kD=k*dk*thickness
        (ok,z0,tfd,tfe)=newton(z0, model, deriv, numRefl)
        if not ok: break
        n[k]=z0
        tfdiff[k]=tfd
        tfexp[k]=tfe
    k_lo =k
# start at kstart and step up
    z0=1.5-0.008J
    for k in range(kstart, M):
        T0=logWfixed[k]
        kD=k*dk*thickness
        (ok,z0,tfd,tfe)=newton(z0, model, deriv, numRefl)
        if not ok: break
        n[k]=z0
        tfdiff[k]=tfd
        tfexp[k]=tfe
    k_hi=k
    return (k_lo, k_hi, n, tfdiff, tfexp)

def convert_to_ri2(kstart, thickness, logWfixed, dk):
    "Fits model to logT to get refractive index"
    global T0, kD
    M=len(logWfixed)
    n=numpy.zeros(M,complex)
# start at good freq and step down
    #z0=2.0-0.005J
    for k in range(kstart, 0, -1):
        T0=logWfixed[k]
        kD=float(k*dk*thickness)
        z0=T0.imag*kD+1.0+(1.0J)*T0.real*kD
        #print z0
        n[k]=z0
    k_lo =k
# start at kstart and step up
    #z0=1.5-0.008J
    for k in range(kstart, M):
        T0=logWfixed[k]
        kD=float(k*dk*thickness)
        z0=T0.imag*kD+1.0+(1.0J)*T0.real*kD
        n[k]=z0
    k_hi=k
    return (k_lo, k_hi, n)

def refine_ri(kstart, thickness, n, model, deriv, W, dk, polarisation_configuration, numRefl_conv):
    "Fits model to raw transmission coefficient to get refractive index"
    "using refractive index values as starting estimates"
    global T0, kD, freqk
    M=len(W)/2
    nrefined=numpy.zeros(M,complex)
    tfdiffconv=numpy.zeros(M,complex)
    tfexpconv=numpy.zeros(M,complex)
# start at good freq and step down
    for k in range(kstart, 0, -1):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        T0=W[k]
        kD=k*dk*thickness
        z0=n[k]
        #(ok,z0,tfdc,tfec)=newtonconv(z0, model, deriv, models.theta_cutoff)
        (ok,z0,tfdc,tfec)=newton_conv(z0, model, deriv, polarisation_configuration, numRefl_conv)
        if not ok: break
        nrefined[k]=z0
        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
    k_lo=k
# start at kstart and step up
    for k in range(kstart, M):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        T0=W[k]
        kD=k*dk*thickness
        z0=n[k]
        #(ok,z0,tfdc,tfec)=newtonconv(z0, model, deriv, models.theta_cutoff)
        (ok,z0,tfdc,tfec)=newton_conv(z0, model, deriv, polarisation_configuration, numRefl_conv)
        if not ok: break
        nrefined[k]=z0
        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
    k_hi=k
    return (k_lo, k_hi, nrefined, tfdiffconv, tfexpconv)


#from scipy.optimize import minimize
#from scipy.optimize import leastsq

def convert_to_theta_leastsq(kstart, thickness, n, model, deriv, W, dk, polarisation_configuration):
    "Fits model to raw transmission coefficient to get refractive index"
    "using refractive index values as starting estimates"
    "Uses least squares method to minimise T0-f(z)=0 *UNFINISHED***"
    global T0, kD, freqk
    M=len(W)/2
    theta0=numpy.zeros(M,complex)
    tfdiffconv=numpy.zeros(M,complex)
    tfexpconv=numpy.zeros(M,complex)  
    
    #Tconvtheta(1.0, nfreq, ops.kD[0], ops.freqk, theta_0, polarisation_configuration) - ops.T0
  
    def residuals(p, k, T0, nfreq, polarisation_configuration):
        #err = model(z, nfreq, polarisation_configuration)
        err = models.Tconvtheta(1.0, nfreq, kD, freqk, p, polarisation_configuration) - T0
        return err
    
# start at good freq and step down 
    for k in range(kstart, 0, -1):
        freqk=k
        nfreq=n[k]
        T0=W[k]
        kD=k*dk*thickness
        p0=0.1
        
        #theta0lsq = minimize(residuals, p0, args=(k, T0, nfreq, polarisation_configuration), method='nelder-mead', options={'xtol': 1e-8, 'disp': False})  
        #theta0[k]=theta0lsq.x
        
        theta0lsq = leastsq(residuals, p0, args=(k, T0, nfreq, polarisation_configuration))
        theta0[k]=theta0lsq[0]
        tfdiffconv[k]=model(theta0[k], nfreq, polarisation_configuration)
        tfexpconv[k]=model(theta0[k], nfreq, polarisation_configuration)+T0
        
    k_lo=k
# start at kstart and step up
    for k in range(kstart, M):
        freqk=k
        nfreq=n[k]
        T0=W[k]
        kD=k*dk*thickness
        p0=0.1
        
        #theta0lsq = minimize(residuals, p0,  args=(k, T0, nfreq, polarisation_configuration), method='nelder-mead', options={'xtol': 1e-8, 'disp': False})  
        #theta0[k]=theta0lsq.x
                
        theta0lsq = leastsq(residuals, p0, args=(k, T0, nfreq, polarisation_configuration))       
        theta0[k]= theta0lsq[0]
        tfdiffconv[k]=model(theta0[k], nfreq, polarisation_configuration)
        tfexpconv[k]=model(theta0[k], nfreq, polarisation_configuration)+T0
        
    k_hi=k
    
    return (k_lo, k_hi, theta0, tfdiffconv, tfexpconv)

def convert_to_theta(kstart, thickness, n, model, deriv, W, dk, polarisation_configuration):
    "Fits model to raw transmission coefficient to get refractive index"
    "using refractive index values as starting estimates"
    global T0, kD, freqk
    M=len(W)/2
    theta0=numpy.zeros(M,complex)
    z0indexfit=numpy.zeros(M,complex)
    tfdiffconv=numpy.zeros(M,complex)
    tfexpconv=numpy.zeros(M,complex)
    numberOfiter=numpy.zeros(M,complex)
    z0=0.1+0.001J
# start at good freq and step down
    for k in range(kstart, 0, -1):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        nfreq=n[k]
        T0=W[k]
        kD=k*dk*thickness        

        (ok,z0,tfdc,tfec,numIter)=newtonconvtheta(z0, model, deriv, nfreq, polarisation_configuration)        
          
        #(ok,z0,tfdc,tfec)=minimiseconvtheta(0.0, 0.6, model, nfreq, polarisation_configuration)
      
        z0index=n[k]
        (ok,z0index,tfdcindex,tfecindex,numIterindex)=newtonconv(z0index, models.convslabfd, models.convslabdashfd, z0, polarisation_configuration) 
           
        if not ok: break
        theta0[k]=z0                 
        z0indexfit[k]=z0index
        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
        numberOfiter[k]=numIter
    k_lo=k
# start at kstart and step up
    z0=theta0[kstart]
    for k in range(kstart+1, M):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)

        freqk=k
        nfreq=n[k]
        T0=W[k]
        kD=k*dk*thickness          
      
        (ok,z0,tfdc,tfec,numIter)=newtonconvtheta(z0, model, deriv, nfreq, polarisation_configuration)
       
        #(ok,z0,tfdc,tfec)=minimiseconvtheta(0.0, 0.6, model, nfreq, polarisation_configuration)
        
        z0index=n[k]
        (ok,z0index,tfdcindex,tfecindex,numIterindex)=newtonconv(z0index, models.convslabfd, models.convslabdashfd, z0, polarisation_configuration) 
        
        if not ok: break
        theta0[k]=z0        
        z0indexfit[k]=z0index
        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
        numberOfiter[k]=numIter
    k_hi=k
    return (k_lo, k_hi, theta0, tfdiffconv, tfexpconv, z0indexfit, numberOfiter)

def convert_to_theta_cone(kstart, thickness, n, model, deriv, W, dk, polarisation_configuration):
    "Fits model to raw transmission coefficient to get refractive index"
    "using refractive index values as starting estimates"
    global T0, kD, freqk
    M=len(W)/2
    theta0=numpy.zeros(M,complex)
    tfdiffconv=numpy.zeros(M,complex)
    tfexpconv=numpy.zeros(M,complex)
    numberOfiter=numpy.zeros(M,complex)
    z0=0.1+0.001J
# start at good freq and step down
    for k in range(kstart, 0, -1):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        nfreq=n[k]
        T0=W[k]
        kD=k*dk*thickness     

        (ok,z0,tfdc,tfec,numIter)=newtonconvtheta(z0, model, deriv, nfreq, polarisation_configuration)        
                           
        if not ok: break
        theta0[k]=z0                 

        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
        numberOfiter[k]=numIter
    k_lo=k
# start at kstart and step up
    z0=theta0[kstart]
    for k in range(kstart+1, M):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)

        freqk=k
        nfreq=n[k]
        T0=W[k]
        kD=k*dk*thickness        
      
        (ok,z0,tfdc,tfec,numIter)=newtonconvtheta(z0, model, deriv, nfreq, polarisation_configuration)
        
        if not ok: break
        theta0[k]=z0        

        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
        numberOfiter[k]=numIter
    k_hi=k
    return (k_lo, k_hi, theta0, tfdiffconv, tfexpconv, numberOfiter)

def fittoreference(z0, model, deriv, k):
    #nreference=readdatafile2('ic2000pnRt.txt', 1.0)
    #nreference=readdatafile2('inRC_00_R_00.txt', 1.0)
    nreference=readdatafile2('i_pt10.0pf8_C_00_R_00.txt', 1.0)
    #nreference=readdatafile2('ec2000pnR.txt', 1.0)
    theta_fdcutoff=0.0001
    (ok,z0,tfdc,tfec)=newtonconv(z0, model, deriv, theta_fdcutoff)
    while z0.real > nreference[k]:
    #while z0.imag > nreference[k]:
        theta_fdcutoff+=0.0001
        #theta_fdcutoff+=0.0001
        (ok,z0,tfdc,tfec)=newtonconv(z0, model, deriv, theta_fdcutoff)
        #print theta_fdcutoff, z0.real, nreference[k]    
    return (ok,z0,tfdc,tfec, theta_fdcutoff)


def refine_ri2(kstart, thickness, n, model, deriv, W, dk):
    "Fits model to raw transmission coefficient to get refractive index"
    "using refractive index values as starting estimates"
    "Beam profiles"
    global T0, kD, freqk
    M=len(W)/2
    fdbeamprofile=numpy.zeros(M,complex)
    nrefined=numpy.zeros(M,complex)
    tfdiffconv=numpy.zeros(M,complex)
    tfexpconv=numpy.zeros(M,complex)
# start at good freq and step down
    for k in range(kstart, 0, -1):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        T0=W[k]
        kD=k*dk*thickness
        z0=n[k]
        (ok,z0,tfdc,tfec, theta_fdcutoff)=fittoreference(z0, model, deriv, k)
        #(ok,z0,tfdc,tfec)=newtonconv(z0, model, deriv)
        if not ok: break
        fdbeamprofile[k]=theta_fdcutoff
        nrefined[k]=z0
        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
        #print k
    k_lo=k
# start at kstart and step up
    for k in range(kstart, M):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        T0=W[k]
        kD=k*dk*thickness
        z0=n[k]
        (ok,z0,tfdc,tfec, theta_fdcutoff)=fittoreference(z0, model, deriv, k)
        #(ok,z0,tfdc,tfec)=newtonconv(z0, model, deriv)
        if not ok: break
        fdbeamprofile[k]=theta_fdcutoff
        nrefined[k]=z0
        tfdiffconv[k]=tfdc
        tfexpconv[k]=tfec
        #print k
    k_hi=k
    return (k_lo, k_hi, nrefined, tfdiffconv, tfexpconv, fdbeamprofile)

def detrend(z, weight, M):
    "Removes linear trend from z; used for estimate of thickness"
    krange = numpy.array(range(0,M))
    sumw=sum(weight[0:M])
    sumwk=sum(weight[0:M]*krange)
    sumwksq=sum(weight[0:M]*krange*krange)
    sumwz=sum(weight[0:M]*z)
    sumwzk=sum(weight[0:M]*krange*z)
    b=(sumwz*sumwksq-sumwzk*sumwk)/(sumw*sumwksq-sumwk*sumwk)
    a=(sumwzk*sumw-sumwz*sumwk)/(sumw*sumwksq-sumwk*sumwk)
    return z-a*krange-b


def PYlevinson(r, y):
# solves sum_{j=0}^{N-1} R_{i-j} x_j = y_i, i=0 ... N-1.
# R_j stored in r[N+j].  r[0] is top right element, r[2N-2] bottom left
    n=len(y)
    g=numpy.zeros(n, float)
    h=numpy.zeros(n, float)
    x=numpy.zeros(n, float)
    x[0]=y[0]/r[n-1]
    g[0]=r[n-2]/r[n-1]
    h[0]=r[n]/r[n-1]
    for m in range(0,n):
        mrange=range(0,m+1)
        m1=m+1
        sxn = -y[m1]
        sd  = -r[n-1]
        for j in mrange:
            sxn += r[n+m-j]*x[j]
            sd  += r[n+m-j]*g[m-j]
        x[m1]=sxn/sd
        for j in mrange: x[j] -= x[m1]*g[m-j]
        if m1==n-1: return x
        sgn = -r[n-m-3]
        shn = -r[n+m+1]
        sgd = -r[n-1]
        for j in mrange:
            sgn += r[n+j-m-2]*g[j]
            shn += r[n+m-j]*h[j]
            sgd += r[n+j-m-2]*h[m-j]
        g[m1]=sgn/sgd
        h[m1]=shn/sd
        k=m
        m2=(m+2)>>1
        pp=g[m1]
        qq=h[m1]
        for j in range(0,m2):
            pt1=g[j]
            pt2=g[k]
            qt1=h[j]
            qt2=h[k]
            g[j]=pt1-pp*qt2
            g[k]=pt2-pp*qt1
            h[j]=qt1-qq*pt2
            h[k]=qt2-qq*pt1
            k-=1

def Clevinson(R,y):
    " wrapper for the Cfunction _levinson"
    n=numpy.size(y)
    x=numpy.zeros(n,float)
    work1=numpy.zeros(n,float)
    work2=numpy.zeros(n,float)
    _levinson.pylev(R, y, x, work1, work2)
    return x

def newtonthinfilm(z0, f, fdash, ns):
    "Solves f(z)=0 by Newton's method, returns tuple (ok, z)"
    maxstep=0.1
    dz=(100.0+100.0J)
    i=0
    z=z0
    ok=FileDialog.TRUE
    while abs(dz)>1e-7:
        dz = -f(z, ns)/fdash(z, ns)
        if abs(dz) > maxstep: dz = dz*maxstep/abs(dz)
        z=z+dz
        i+=1
        if i>100:
            ok=FileDialog.FALSE
            break
    return (ok, z, f(z, ns),f(z, ns)+T0)

def newton(z0, f, fdash, numRefl):
    "Solves f(z)=0 by Newton's method, returns tuple (ok, z)"
    maxstep=0.1
    dz=(100.0+100.0J)
    i=0
    z=z0
    ok=FileDialog.TRUE
    while abs(dz)>1e-7:
        dz = -f(z, numRefl)/fdash(z, numRefl)
        if abs(dz) > maxstep: dz = dz*maxstep/abs(dz)
        z=z+dz
        i+=1
        if i>100:
            ok=FileDialog.FALSE
            break
    return (ok, z, f(z, numRefl),f(z, numRefl)+T0)

def newton_conv(z0, f, fdash, polarisation_configuration, numRefl_conv):
    "Solves f(z)=0 by Newton's method, returns tuple (ok, z)"
    maxstep=0.1
    dz=(100.0+100.0J)
    i=0
    z=z0
    ok=FileDialog.TRUE
    while abs(dz)>1e-7:
        dz = -f(z, polarisation_configuration, numRefl_conv)/fdash(z, polarisation_configuration, numRefl_conv)
        if abs(dz) > maxstep: dz = dz*maxstep/abs(dz)
        z=z+dz
        i+=1
        if i>100:
            ok=FileDialog.FALSE
            break
    return (ok, z, f(z, polarisation_configuration, numRefl_conv),f(z, polarisation_configuration, numRefl_conv)+T0)

def newtonconv(z0, f, fdash, theta_fdcutoff, polarisation_configuration):
    "Solves f(z)=0 by Newton's method, returns tuple (ok, z)"
    maxstep=0.1
    dz=(100.0+100.0J)
    i=0
    z=z0
    ok=FileDialog.TRUE
    while abs(dz)>1e-7:
        dz = -f(z,theta_fdcutoff, polarisation_configuration)/fdash(z,theta_fdcutoff, polarisation_configuration)
        if abs(dz) > maxstep: dz = dz*maxstep/abs(dz)
        z=z+dz
        i+=1
        if i>100:
            ok=FileDialog.FALSE
            break
    return (ok, z, f(z,theta_fdcutoff, polarisation_configuration),f(z,theta_fdcutoff, polarisation_configuration)+T0,i)

def newtonconvtheta(z0, f, fdash, nfreq, polarisation_configuration):
    "Solves f(z)=0 by Newton's method, returns tuple (ok, z)"
    maxstep=0.1
    dz=(100.0+100.0J)
    i=0
    z=z0
    ok=FileDialog.TRUE
    while abs(dz)>1e-7:
        dz = -f(z,nfreq,  polarisation_configuration)/fdash(z,nfreq, polarisation_configuration)
        if abs(dz) > maxstep: dz = dz*maxstep/abs(dz)
        z=z+dz
        i+=1
        if i>100:
            ok=FileDialog.FALSE
            break
    return (ok, z, f(z,nfreq, polarisation_configuration),f(z,nfreq, polarisation_configuration)+T0, i)

def minimiseconvtheta(a, b, f, nfreq, polarisation_configuration):
    "Finds the minimum of f(z)=0"
    
    funct=numpy.zeros(100,complex)
    base=numpy.zeros(100,complex)
    k=0
    for i in numpy.linspace(0,1.0,100):
        funct[k]=f(i, nfreq, polarisation_configuration)
        base[k]=i
        k+=1
     
    #pylab.figure()
    pylab.plot(base[1:], funct[1:].real)
    #pylab.plot(base[1:], funct[1:].imag)
    #pylab.show()
     
# This function takes the minimum value of the function f(x) which is 
#  assumed to lie between a and b
    r=0.618034
    delta=1e-4
    x=a+r*r*(b-a); y=a+r*(b-a)
    fx=f(x, nfreq, polarisation_configuration); fy=f(y, nfreq, polarisation_configuration)
    ok=FileDialog.TRUE

    while True:
        if fy<fx:

            a=x; x=y; y=a+r*(b-a); fx=fy; fy=f(y, nfreq, polarisation_configuration)

     
        else:

            b=y; y=x; x=a+r*r*(b-a); fy=fx; fx=f(x, nfreq, polarisation_configuration)

           
        if abs(y-x)<delta*abs(x): break
    

    return (ok, x, f(x,nfreq, polarisation_configuration),f(x,nfreq, polarisation_configuration)+T0)


