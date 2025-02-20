from numpy import *
import etc
from numpy.fft import fft, ifft
from physics import graph
#from cmath import *      # note this prevents Ufuncs existing
#
global T0, kD, freqk
T0=0.0; kD=[0.0]

def read_raw_datafile(cols, filename, scale):
    "reads datafile, ignoring lines beginning with #, returns arrays"
    "containing reading and error, computes sampling interval"
    "if <3 columns, zero sampling interval and null arrays are returned"
    z1=[]; z2=[]
    fileobj=open(filename, 'r')
    for line in fileobj.readlines():
# ignore lines beginning with #
        if line[0]=='#': continue
        if cols==1: z1.append(scale*float(line.split(',')[0]))
        if cols==2:
#            print(line)
            print(line.split('\t'))
            z1.append(scale*float(line.split(',')[1]))
        if cols==3:
            z1.append(scale*float(line.split(',')[1]))
            z2.append(scale*float(line.split(',')[2]))
# try to get sampling interval from last time in file
    sampling_interval=0.0
    if cols>=2: sampling_interval=float(line.split(',')[0])/(len(z1)-1)

    fileobj.close()
    return (sampling_interval, array(z1), array(z2))

def chop_shift_data(z, start, offset, L):
    "Ignores first beginning, embeds in length L, truncating if necessary"
    z1=zeros(L, float)
    l=len(z)
    if l-start+offset<=L:
# needs zero padding
        z1[offset:l-start+offset]=z[start:l]
    else:
# needs truncating
        z1[offset:L]=z[start:start+L-offset]
    return z1

def s(t, N, P, Q):
    "interpolation function"
    c=cos(pi/(Q+2))
    eps=1.12345e-10
    f1=sin((pi*(P-Q)*(t+eps))/N)/sin((pi*(t+eps))/N)
    f2=cos((pi*(Q+2)*(t+eps))/N)/(cos((2*pi*(t+eps))/N) - c)
    return (1.0-c)*f1*f2*cos(pi*t/N)/N

def interpolate_new(x, dt_in, dt_out):
    N=len(x)
    Q=int(0.05*N)
    M=N-2*Q
    eps=(dt_out - dt_in)/dt_in
# band limit to -M/2 to M/2 with smooth transition X[+-M/2]=0
    X=fft(x)
    X[M/2+1:N-M/2]=0.0
    for k in range(-M/2,-M/2+Q):   X[k]=X[k]*0.5*(1.0+cos((pi*(k+M/2-Q))/Q))
    for k in range( M/2-Q, M/2+1): X[k]=X[k]*0.5*(1.0+cos((pi*(k-M/2+Q))/Q))
    x=ifft(X).real
# interpolate
    P=N
    Q=Q-1
    R=10
    z=zeros(N, float)
    for m in range(0,N):
        a=0.0
        T=m*dt_out/dt_in
        n0=int(T)
        for n in range(n0-R, n0+R):
            a=a+x[n%N]*s(T - n, N, P, Q)

        z[m]=a
    return z

def interpolate(x, dt_in, dt_out):
    N=len(x)
    Q=int(0.05*N)
    M=N-2*Q
    eps=(dt_out - dt_in)/dt_in
# band limit to -M/2 to M/2 with smooth transition X[+-M/2]=0
    X=fft(x)
    X[M/2+1:N-M/2]=0.0
    for k in range(-M/2,-M/2+Q):   X[k]=X[k]*0.5*(1.0+cos((pi*(k+M/2-Q))/Q))
    for k in range( M/2-Q, M/2+1): X[k]=X[k]*0.5*(1.0+cos((pi*(k-M/2+Q))/Q))
    x=ifft(X).real
# interpolate
    P=N
    Q=Q-1
    R=10
    z=zeros(N, float)
    for m in range(0,N):
        a=0.0
        for r in range(-R, R):
            a=a+x[(m+r)%N]*s(r-m*eps, N, P, Q)

        z[m]=a
    return z

def decimate(x, L, newL):
    "converts data array x of L points to zero-padded array"
    "newL points to reduce oversampling, both L and newL must be powers of 2"
    X=fft(x[0:L])
    Z=zeros(newL, complex)
    LL=newL/2
    Z[0:LL]=X[0:LL]; Z[-LL:L]=X[-LL:L]
    x=ifft(Z).real
    return x

def ignore_data(dy, ignore):
    "marks sections of data to be ignored by 'infinite' std dev"
    dy1=array(dy)
    for (start, finish) in ignore: dy1[start:finish]=1e30
    return dy1

def savedatafile(cols, deltat, z1, z2, filename, scale):
    "saves data files.  Make z2=z1 if not used"
    L=len(z1)
    fileobj=open(filename, 'w')
    for i in range(0,L):
        if cols==1: fileobj.write(str(z1*scale) + '\n')
        if cols==2: fileobj.write(str(i*deltat)+'  ' + str(z1[i]*scale) + '\n')
        if cols==3:
            fileobj.write(str(i*deltat)    + '  ' +   \
                          str(z1[i]*scale) + '  ' +   \
                          str(z2[i]*scale) + '\n')
    fileobj.close()

def applyEXX(z):
    "used in conjgrad method of getting impulse response"
    M=len(z); N=2*M
    w=zeros(N, float)
    w[0:M]=z[0:M]
    W=fft(w)/N
    W=N*W*globXX
    w=N*ifft(W).real
    W=fft(globee*w)/N
    W=W*globXXconj
    return (N*N*ifft(W).real)[0:M]

def chilim(w):
    "returns True if chisq < data length, used to stop fitting"
    M=len(w); N=2*M
    z=zeros(N, float)
    z[0:M]=w[0:M]
    Z=fft(z)/N
    W=N*Z*globXX
    ypred=N*ifft(W).real
    return sum(globee*(globyy - ypred)**2)<M

def impulse_fit(XX, XXconj, yy, ee, matrix, rhs, w_init):
    "The only reason for this is to use a better name and keep etc hidden"
    global globXX, globXXconj, globyy, globee
    globXX=XX; globXXconj=XXconj; globyy=yy; globee=ee
    return etc.conjgrad(matrix, rhs, w_init, chilim, 200, False)[0]

def truncate(z, length, fall):
    "Removes end of data with smooth cut-off"
    N=len(z)
# find peak
    peak=-1e10; place=0
    for t in range(0,N):
        if abs(z[t])>peak: peak=abs(z[t]); place=t

# set weight function
    weight=time_window(length+place, fall, N)
    return z*weight

def time_window(where, scale, N):
    "model for removing data with long time delays"
    weight=zeros(N,float)
    for t in range(0,N):
        weight[t]=1.0/(1.0 + exp((t-where)/float(scale)))
    return weight

def noise_model_time(z):
    "Totally empirical model to get noise estimate"
    N=len(z)
    Z=fft(abs(z))/N
    K=100
    for i in range(1,K):
        f=exp(-0.001*i*i)
        Z[i]=Z[i]*f
        Z[N-i]=Z[N-i]*f
    Z[K:N-K]=0.0
    z=N*ifft(Z).real
    return 5e-4 + 0.02*z*z

def good_signal(N, kstart, goodrange):
    "Model of range where signal/noise is good, used for fixing phase"
    weight=zeros(N, float)
    M=N/2
    for k in range(10,M+1):
        weight[k]=exp(-((k-kstart)/float(goodrange))**2)
        weight[N-k]=weight[k]
    weight[0]=0.0;
    return weight

def fixphase(logW, kstart, weight):
    "unwraps phase, starting at freq where S/N is good and working outwards"
    M=len(logW)
    phase=array(logW.imag)
# start at freq kstart and step down
    lastph=phase[kstart]
    n=0
    for k in range(kstart, -1, -1):
        if phase[k]-lastph >  1.4*pi: n=n-1
        if phase[k]-lastph < -1.4*pi: n=n+1
        lastph=phase[k]
        phase[k]=phase[k] + 2*n*pi
# start at freq kstart and step up
    lastph=phase[kstart]
    n=0
    for k in range(kstart, M, 1):
        if phase[k]-lastph >  1.4*pi: n=n-1
        if phase[k]-lastph < -1.4*pi: n=n+1
        lastph=phase[k]
        phase[k]=phase[k] + 2*n*pi
# extrapolate to zero frequency
    krange = array(range(0,M))
    sumw=sum(weight[0:M])
    sumwk=sum(weight[0:M]*krange)
    sumwksq=sum(weight[0:M]*krange*krange)
    sumwph=sum(weight[0:M]*phase)
    sumwphk=sum(weight[0:M]*krange*phase)
    b=(sumwph*sumwksq-sumwphk*sumwk)/(sumw*sumwksq-sumwk*sumwk)
    a=(sumwphk*sumw-sumwph*sumwk)/(sumw*sumwksq-sumwk*sumwk)
    p=round(b/(2*pi))
    phase=phase-2*p*pi
    return logW.real + (1J)*phase

def convert_to_ri(kstart, thickness, model, deriv, logWfixed, dk):
    "Fits model to logT to get refractive index"
    global T0, kD
    M=len(logWfixed)
    n=zeros(M,complex)
# start at good freq and step down
    z0=1.5-0.008J
    for k in range(kstart, 0, -1):
        T0=logWfixed[k]
        kD=k*dk*thickness
        (ok,z0)=newton(z0, model, deriv)
        if not ok: break
        n[k]=z0
    k_lo =k
# start at kstart and step up
    z0=1.5-0.008J
    for k in range(kstart, M):
        T0=logWfixed[k]
        kD=k*dk*thickness
        (ok,z0)=newton(z0, model, deriv)
        if not ok: break
        n[k]=z0
    k_hi=k
    return (k_lo, k_hi, n)

def refine_ri(kstart, thickness, n, model, deriv, W, dk):
    "Fits model to raw transmission coefficient to get refractive index"
    "using refractive index values as starting estimates"
    global T0, kD, freqk
    M=len(W)/2
    nrefined=zeros(M,complex)
# start at good freq and step down
    for k in range(kstart, 0, -1):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        T0=W[k]
        kD=k*dk*thickness
        z0=n[k]
        (ok,z0)=newton(z0, model, deriv)
        if not ok: break
        nrefined[k]=z0
    k_lo=k
# start at kstart and step up
    for k in range(kstart, M):
#        T0=complex((log(Y[k]/X[k])).real, logWfixed[k].imag)
        freqk=k
        T0=W[k]
        kD=k*dk*thickness
        z0=n[k]
        (ok,z0)=newton(z0, model, deriv)
        if not ok: break
        nrefined[k]=z0
    k_hi=k
    return (k_lo, k_hi, nrefined)

def detrend(z, weight, M):
    "Removes linear trend from z; used for estimate of thickness"
    krange = array(range(0,M))
    sumw=sum(weight[0:M])
    sumwk=sum(weight[0:M]*krange)
    sumwksq=sum(weight[0:M]*krange*krange)
    sumwz=sum(weight[0:M]*z)
    sumwzk=sum(weight[0:M]*krange*z)
    b=(sumwz*sumwksq-sumwzk*sumwk)/(sumw*sumwksq-sumwk*sumwk)
    a=(sumwzk*sumw-sumwz*sumwk)/(sumw*sumwksq-sumwk*sumwk)
    return z-a*krange-b

def newton(z0, f, fdash):
    "Solves f(z)=0 by Newton's method, returns tuple (ok, z)"
    maxstep=0.1
    dz=(1.0+1.0J)
    i=0
    z=z0
    ok=True
    while abs(dz)>1e-7:
        dz = -f(z)/fdash(z)
        if abs(dz) > maxstep: dz = dz*maxstep/abs(dz)
        z=z+dz
        i+=1
        if i>100:
            ok=False
            break
    return (ok, z)

