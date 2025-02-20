import numpy
# import FileDialog
import pylab

global T0, kD, freqk
T0=0.0; kD=[0.0]
    

def logT3(n0, n1, kD):
    return numpy.log(4.0*n0*n1/(n0+n1)**2) - (1J)*(n1-1.0)*kD

def logT3dash(n0, n1, kD):
    return 1.0/n1 - 2.0/(n0+n1) - (1J)*kD

def logT4(n0, n1, kD):
    return numpy.log(4.0*n0*n1/(n0+n1)**2) - (1J)*(n1-1.0)*kD

def logT4dash(n0, n1, kD):
    return 1.0/n1 - 2.0/(n0+n1) - (1J)*kD

def newfixphase(W, kstart, weight, M):
    "unwraps phase using phase differences"
    logW = numpy.log(W[0:M])
    M=len(logW)
    phase=numpy.zeros(M, float)
        
    phase[kstart] = logW[kstart].imag
    print kstart, M    
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
    #phase=phase-b#2*p*numpy.pi
    phase = phase-2*p*numpy.pi        
    
    return logW.real + (1J)*(phase )




def convert_to_ri(kstart, thickness, model, deriv, logWfixed, dk):
    global T0, kD
    M=len(logWfixed)
    n=numpy.zeros(M,complex)
    tfdiff=numpy.zeros(M,complex)
    tfexp=numpy.zeros(M,complex)

    z0=2.0-0.005J
    
    for k in range(kstart, 0, -1):
        T0=logWfixed[k]
        kD=k*dk*thickness
        (ok,z0,tfd,tfe)=newton(z0, model, deriv)
        if not ok: break
        n[k]=z0
        tfdiff[k]=tfd
        tfexp[k]=tfe
    k_lo =k

    z0=2.0-0.005J
    for k in range(kstart, M):
        T0=logWfixed[k]
        kD=k*dk*thickness
        (ok,z0,tfd,tfe)=newton(z0, model, deriv)
        if not ok: break
        n[k]=z0
        tfdiff[k]=tfd
        tfexp[k]=tfe
    k_hi=k
    return (k_lo, k_hi, n, tfdiff, tfexp)

def newton(z0, f, fdash):
    maxstep=0.1
    dz=(100.0+100.0J)
    i=0
    z=z0
    # ok=FileDialog.TRUE
    while abs(dz)>1e-13:
        dz = -f(z)/fdash(z)
        if abs(dz) > maxstep: dz = dz*maxstep/abs(dz)
        z=z+dz
        i+=1
        if i>100:
            # ok=FileDialog.FALSE
            break
        ok=1
    return (ok, z, f(z),f(z)+T0)

