# Functions for models included here.
#
# T1(n0, n1, kD)         Transmission coefficient of single layer
# R1(n0, n1, kD)         Reflection coefficient of single layer
# rts(n, kD)             Tuple of r, t and t**2 - r**2 for single layer"
# RTm(m, layers)         Transmission coefficient of m layers
#
# T1dash(n0, n1, kD)     Derivative of T1 wrt n1
# R1dash(n0, n1, kD)     Derivative of R1 wrt n1
# rtsdash(n, kD)         Derivatives of r,t,s wrt n
# RTmdash(m, layers, q)  Derivative of T for m layers wrt index of layer q
#
# logT1(n0, n1, kD)              log of transmission coefficient for 1 layer
# logT2(n0, n1, ,n2, kD1, kD2)   log of transmission coefficient for 2 layers
# logT3(n0,n1,n2,n3,kD1,kD2,kD3) log of transmission coefficient for 3 layers
# logTm(m, layers)               log of transmission coefficient for m layers
#
# logT1dash(n0, n1, kD)           derivative of logT1 wrt n1
# logT2dash(n0, n1, n2, kD1, kD2) derivative of logT2 wrt n2
# logT3dash(n0,n1,n2,kD1,kD2,kD3) derivative of logT3 wrt n2
#
# Note: Reflection coefficient for more than one layer are for waves incident
#       from the EXIT end.
#       The phase of all transmission coefficients is relative to the same
#       total thickness of material of index n0.
import pylab
import numpy
from numpy.fft import fft, ifft
#from cmath import *      # note this prevents Ufuncs existing
import ops
import random


#width_file="w5m5552m2025.txt"
#width=ops.readdatafilewidth(width_file, 1.0)
#width2=numpy.zeros(len(width),float)

#for i in range(0,len(width)):
#    inf=width[2000]
#    nan=width[4]
#    if width[i] == nan:
#        width2[i]=0.0 
#    elif width[i] == inf:
#        width2[i]=0.0 
#    else:
#        width2[i]=width[i]

#
global T0, kD, kDs
T0=0.0; kD=0.0; kDs=0.0

    # Single slab of refractive index n in free space
def rawslab1(n):
    thickscale=1.0
    return T1(1.0, n, thickscale*ops.kD[0]) - ops.T0

def rawslab1dash(n):
    thickscale=1.0
    return T1dash(1.0, n, thickscale*ops.kD[0])

def slab1(n, numRefl):
    thickscale=1.0
    return logT1(1.0, n, thickscale*ops.kD[0], numRefl) - ops.T0

def slab1dash(n, numRefl):
    thickscale=1.0
    return logT1dash(1.0, n, thickscale*ops.kD[0], numRefl)

def thinfilm1(n, ns):
    thickscale=1.0
    return logthinfilmT1(1.0, n, thickscale*ops.kD[0], ns, kDs) - ops.T0

def thinfilm1dash(n, ns):
    thickscale=1.0
    return logthinfilmT1dash(1.0, n, thickscale*ops.kD[0], ns, kDs)

def logthinfilmT1(n0, n1, kD, ns, kDs):
    "logT for single layer with NO REFLECTION"
    return numpy.log(2.0*n1*(n0+ns)/((n0+n1)*(n1+ns))) - (1J)*(n1-n0)*kD
    "logT for single layer with inf REFLECTION"
    #return numpy.log(2.0*n1*(n0+ns)/((n0+n1)*(n1+ns))) - (1J)*(n1-n0)*kD - numpy.log(1.0-thinfilmE1(n0, n1, kD, ns, kDs)) - numpy.log(1.0-thinfilmE2(n0, n1, kD, ns, kDs)) + numpy.log(1.0-thinfilmE3(n0, n1, kD, ns, kDs))
   
def thinfilmE1(n0, n1, kD, ns, kDs):
    return ((n1-ns)/(n1+ns))*((n1-n0)/(n1+n0))*numpy.exp((-2.0J)*n1*kD)

def thinfilmE2(n0, n1, kD, ns, kDs):
    return ((ns-n0)/(ns+n0))*((ns-n1)/(ns+n1))*numpy.exp((-2.0J)*ns*kDs)

def thinfilmE3(n0, n1, kD, ns, kDs):
    return ((ns-n0)/(ns+n0))**2*numpy.exp((-2.0J)*ns*kDs)

def logthinfilmT1dash(n0, n1, kD, ns, kDs):
    "Derivative of logT for single layer with NO REFLECTION"
    return 1.0/n1 - 1.0/(n0+n1) - 1.0/(n1+ns) - (1J)*kD
    "Derivative of logT for single layer with inf REFLECTION"
    #return 1.0/n1 - 1.0/(n0+n1) - 1.0/(n1+ns) - (1J)*kD + thinfilmE1dash(n0, n1, kD, ns, kDs)/thinfilmE1(n0, n1, kD, ns, kDs) + thinfilmE2dash(n0, n1, kD, ns, kDs)/thinfilmE2(n0, n1, kD, ns, kDs) - thinfilmE3dash(n0, n1, kD, ns, kDs)/thinfilmE3(n0, n1, kD, ns, kDs)
    
def thinfilmE1dash(n0, n1, kD, ns, kDs):
    return thinfilmE1(n0, n1, kD, ns, kDs)*(1.0/(n1-ns) - 1.0/(n1+ns) + 1.0/(n1-n0) - 1.0/(n1+n0) + (-2.0J)*kD)
    
def thinfilmE2dash(n0, n1, kD, ns, kDs):
    return thinfilmE2(n0, n1, kD, ns, kDs)*(1.0/(ns-n0) - 1.0/(ns+n0) + 1.0/(ns-n1) - 1.0/(ns+n1) + (-2.0J)*kDs)

def thinfilmE3dash(n0, n1, kD, ns,  kDs):
    return thinfilmE3(n0, n1, kD, ns, kDs)*(2.0/(ns-n0) - 2.0/(ns+n0) + (-2.0J)*kDs)

def slab1para(n):
    thickscale=1.0
    return logT1para(1.0, n, thickscale*ops.kD[0], sampleangle) - ops.T0

def slab1dashpara(n):
    thickscale=1.0
    return logT1dashpara(1.0, n, thickscale*ops.kD[0], sampleangle)

def slab1perp(n):
    thickscale=1.0
    return logT1perp(1.0, n, thickscale*ops.kD[0], sampleangle) - ops.T0

def slab1dashperp(n):
    thickscale=1.0
    return logT1dashperp(1.0, n, thickscale*ops.kD[0], sampleangle)

# Three layer sandwich, fitting middle layer
def GBsand(nfill):
    thickscale=1.0
    return logT3(1.0, nbread, nfill, nbread, \
                 ops.kD[0], thickscale*ops.kD[1], ops.kD[2]) - ops.T0

def GBsandddkash(nfill):
    thickscale=1.0
    return logT3dash(1.0, nbread, nfill, nbread, \
                     ops.kD[0], thickscale*ops.kD[1], ops.kD[2])

# Filled roll, fitting outer layers
def filled_roll(nbread):
    thickscale=1.0
    return logT3(1.0, nbread, nfill, nbread, \
                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]) - ops.T0

def filled_rolldash(nbread):
    thickscale=1.0
    dn=0.01
    return (logT3(1.0, nbread+dn, nfill, nbread+dn, \
                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]) - \
            logT3(1.0, nbread-dn, nfill, nbread-dn, \
                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]))/(2.0*dn)

# Single slab of refractive index n observed using converging beam
# note this is T not logT and is used only for refinement
def convslab(n, polarisation_configuration, numRefl_conv):
    return Tconv(1.0, n, ops.kD[0], ops.freqk, polarisation_configuration, numRefl_conv) - ops.T0

def convslabdash(n, polarisation_configuration, numRefl_conv):
    dn=0.01 + 0.001J
    return (Tconv(1.0, n+dn, ops.kD[0], ops.freqk, polarisation_configuration, numRefl_conv) - \
            Tconv(1.0, n-dn, ops.kD[0], ops.freqk, polarisation_configuration, numRefl_conv))/(2.0*dn) 
            
def convslabtheta(theta_0,nfreq, polarisation_configuration):
    return Tconvtheta(1.0, nfreq, ops.kD[0], ops.freqk, theta_0, polarisation_configuration)- ops.T0

def convslabdashtheta(theta_0,nfreq, polarisation_configuration):
    dtheta_0=0.01 + 0.001J
    return (Tconvtheta(1.0, nfreq, ops.kD[0], ops.freqk, theta_0+dtheta_0, polarisation_configuration) - \
            Tconvtheta(1.0, nfreq, ops.kD[0], ops.freqk, theta_0-dtheta_0, polarisation_configuration))/(2.0*dtheta_0)
            
def convslabthetacone(theta_0,nfreq, polarisation_configuration):
    return Tconvthetacone(1.0, nfreq, ops.kD[0], ops.freqk, theta_0, polarisation_configuration) - ops.T0

def convslabdashthetacone(theta_0,nfreq, polarisation_configuration):
    dtheta_0=0.01 + 0.001J
    return (Tconvthetacone(1.0, nfreq, ops.kD[0], ops.freqk, theta_0+dtheta_0, polarisation_configuration) - \
            Tconvthetacone(1.0, nfreq, ops.kD[0], ops.freqk, theta_0-dtheta_0, polarisation_configuration))/(2.0*dtheta_0)
            
def convslabdonut(n, polarisation_configuration):
    return Tconvdonut(1.0, n, ops.kD[0], ops.freqk, polarisation_configuration) - ops.T0

def convslabdashdonut(n, polarisation_configuration):
    dn=0.01 + 0.001J
    return (Tconvdonut(1.0, n+dn, ops.kD[0], ops.freqk, polarisation_configuration) - \
            Tconvdonut(1.0, n-dn, ops.kD[0], ops.freqk, polarisation_configuration))/(2.0*dn) 

def convslabfd(n,theta_fdcutoff, polarisation_configuration):
    return Tconvfd(1.0, n, ops.kD[0], ops.freqk, theta_fdcutoff,polarisation_configuration) - ops.T0

def convslabdashfd(n,theta_fdcutoff, polarisation_configuration):
    dn=0.01 + 0.001J
    return (Tconvfd(1.0, n+dn, ops.kD[0], ops.freqk, theta_fdcutoff, polarisation_configuration) - \
            Tconvfd(1.0, n-dn, ops.kD[0], ops.freqk, theta_fdcutoff, polarisation_configuration))/(2.0*dn)
            
def logT1para(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logT for single layer with inf reflections for electric field parallel to the plane of incidence"
    return numpy.log(4.0*n0*n1*costheta*costheta_t/(n0*costheta_t+n1*costheta)**2)-(1J)*(n1*costheta_t-n0*costheta)*kD - numpy.log(1.0-E1para(n0,n1,kD,theta))

def E1para(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return ((n0*costheta_t-n1*costheta)/(n1*costheta+n0*costheta_t))**2*numpy.exp((-2.0J)*n1*kD/costheta_t)

def E1dashpara(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return E1para(n0, n1, kD, theta)*(((-4.0*n0*costheta_t*costheta)/((n0*costheta_t)**2-(n1*costheta)**2))-2.0J*kD/costheta_t)
    
def logT1dashpara(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logTdash for single layer with inf reflections for electric field parallel to the plane of incidence"
    return 1.0/n1-2.0*costheta/(n0*costheta_t+n1*costheta)**2-(1.0J)*kD*costheta_t + E1dashpara(n0, n1, kD, theta)/(1.0 - E1para(n0, n1, kD, theta))

def logT1perp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logT for single layer with inf reflections for electric field perpendicular to the plane of incidence"
    return numpy.log(4.0*n0*n1*costheta*costheta_t/(n0*costheta+n1*costheta_t)**2)-(1J)*(n1*costheta_t-n0*costheta)*kD - numpy.log(1.0-E1perp(n0,n1,kD,theta))

def E1perp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return ((n1*costheta_t-n0*costheta)/(n1*costheta_t+n0*costheta))**2*numpy.exp((-2.0J)*n1*kD/costheta_t)

def E1dashperp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return E1perp(n0, n1, kD, theta)*(((4.0*n0*costheta_t*costheta)/((n1*costheta_t)**2-(n0*costheta)**2))-2.0J*kD/costheta_t)
    
def logT1dashperp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logTdash for single layer with inf reflections for electric field perpendicular to the plane of incidence"
    return 1.0/n1-2.0*costheta_t/(n0*costheta+n1*costheta_t)**2-(1.0J)*kD*costheta_t + E1dashperp(n0, n1, kD, theta)/(1.0 - E1perp(n0, n1, kD, theta))
    
def logT1para(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logT for single layer with inf reflections for electric field parallel to the plane of incidence"
    return numpy.log(4.0*n0*n1*costheta*costheta_t/(n0*costheta_t+n1*costheta)**2)-(1J)*(n1*costheta_t-n0*costheta)*kD - numpy.log(1.0-E1para(n0,n1,kD,theta))

def E1para(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return ((n0*costheta_t-n1*costheta)/(n1*costheta+n0*costheta_t))**2*numpy.exp((-2.0J)*n1*kD/costheta_t)

def E1dashpara(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return E1para(n0, n1, kD, theta)*(((-4.0*n0*costheta_t*costheta)/((n0*costheta_t)**2-(n1*costheta)**2))-2.0J*kD/costheta_t)
    
def logT1dashpara(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logTdash for single layer with inf reflections for electric field parallel to the plane of incidence"
    return 1.0/n1-2.0*costheta/(n0*costheta_t+n1*costheta)**2-(1.0J)*kD*costheta_t + E1dashpara(n0, n1, kD, theta)/(1.0 - E1para(n0, n1, kD, theta))

def logT1perp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logT for single layer with inf reflections for electric field perpendicular to the plane of incidence"
    return numpy.log(4.0*n0*n1*costheta*costheta_t/(n0*costheta+n1*costheta_t)**2)-(1J)*(n1*costheta_t-n0*costheta)*kD - numpy.log(1.0-E1perp(n0,n1,kD,theta))

def E1perp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return ((n1*costheta_t-n0*costheta)/(n1*costheta_t+n0*costheta))**2*numpy.exp((-2.0J)*n1*kD/costheta_t)

def E1dashperp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    return E1perp(n0, n1, kD, theta)*(((4.0*n0*costheta_t*costheta)/((n1*costheta_t)**2-(n0*costheta)**2))-2.0J*kD/costheta_t)
    
def logT1dashperp(n0, n1, kD, theta):
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    "logTdash for single layer with inf reflections for electric field perpendicular to the plane of incidence"
    return 1.0/n1-2.0*costheta_t/(n0*costheta+n1*costheta_t)**2-(1.0J)*kD*costheta_t + E1dashperp(n0, n1, kD, theta)/(1.0 - E1perp(n0, n1, kD, theta))

def T1(n0, n1, kD):
    "Transmission coefficient of single layer"
    a=4.0*n0*n1/(n0+n1)**2
    b=((n0-n1)/(n0+n1))**2
    return a*numpy.exp(-(1.0J)*(n1-n0)*kD)/(1.0 - b*numpy.exp(-(2.0J)*n1*kD))

def R1(n0, n1, kD):
    "Reflection coefficient of single layer"
    return ((n0-n1)/(n0+n1))*(1-numpy.exp((-2.0J)*n1*kD))/(1-E1(n0, n1, kD))

def T1dash(n0, n1, kD):
    "Derivative of T1"
    return T1(n0, n1, kD)*logT1dash(n0, n1, kD)

def R1dash(n0, n1, kD):
    "Derivative of R1"
    a=-1.0/(n0-n1) - 1.0/(n0+n1)
    b=(2.0J)*kD/(numpy.exp((2.0J)*n1*kD) - 1.0)
    return R1(n0, n1, kD)*(a + b + E1dash(n0, n1, kD)/(1.0 - E1(n0, n1,kD)))

def logT1(n0, n1, kD, numRefl):
    if numRefl == 'none':
        "logT for single layer with NO REFLECTION"
        return numpy.log(4.0*n0*n1/(n0+n1)**2) - (1J)*(n1-n0)*kD

    if  numRefl == 'inf':
        "logT for single layer with inf REFLECTION"
        return numpy.log(4.0*n0*n1/(n0+n1)**2) - (1J)*(n1-n0)*kD - numpy.log(1.0-E1(n0, n1, kD))

    if  numRefl == 'one':
        "logT for single layer with ONE REFLECTION"
        return numpy.log(4.0*n0*n1/(n0+n1)**2) - (1J)*(n1-n0)*kD + numpy.log(1.0+E1(n0, n1, kD)) 

def E1(n0, n1, kD):
    return ((n1-n0)/(n1+n0))**2*numpy.exp((-2.0J)*n1*kD)

def logT1dash(n0, n1, kD, numRefl):
    
    if  numRefl == 'none':
        "Derivative of logT for single layer with NO REFLECTION"
        return 1.0/n1 - 2.0/(n0+n1) - (1J)*kD

    if  numRefl == 'inf':
        "Derivative of logT for single layer with inf REFLECTION"
        return 1.0/n1 - 2.0/(n0+n1) - (1J)*kD + E1dash(n0, n1, kD)/(1.0 - E1(n0, n1,kD))
    
    if  numRefl == 'one':
        "Derivative of logT for single layer with ONE REFLECTION"
        return 1.0/n1 - 2.0/(n0+n1) - (1J)*kD + E1dash(n0, n1, kD)/(1.0 + E1(n0, n1,kD))

def E1dash(n0, n1, kD):
    return E1(n0, n1, kD)*(2.0/(n1-n0) - 2.0/(n1+n0) + (-2.0J)*kD)

def logT2(n0, n1, n2, kD1, kD2):
    "logT for two layers"
    return logT1(n0, n1, kD1) + logT1(n0, n2, kD2) - \
           numpy.log(1 - R1(n0, n1, kD1)*R1(n0, n2, kD2))

def logT2dash(n0, n1, n2, kD1, kD2):
    "Derivative of logT2 wrt n2"
    return logT1dash(n0, n2, kD2) + \
           R1(n0,n1,kD1)*R1dash(n0,n2,kD2)/(1 - R1(n0,n1,kD1)*R1(n0,n2,kD2))

def logT3(n0, n1, n2, n3, kD1, kD2, kD3):
    "logT for three layers"
    r0=R1(n0, n1, kD1)
    r2=R1(n0, n3, kD3)
    (r1,t1,s1)=rts(n0, n2, kD2)
    return logT1(n0, n1, kD1) + logT1(n0, n2, kD2) + logT1(n0, n3, kD3) -  \
           numpy.log(1 - r1*r0 - r2*r1 - s1*r2*r0)

def logT3dash(n0,n1,n2,n3,kD1,kD2,kD3):
    "derivative of logT3 wrt n2, obtained by differencing"
    dn=0.01
    return (logT3(n0, n1, n2+dn, n3, kD1, kD2, kD3) - \
            logT3(n0, n1, n2-dn, n3, kD1, kD2, kD3))/(2*dn)

def RTm(m, n0, layers):
    "Reflection and Transmission coefficients for m layers"
    U=0.0; V=1.0; T=1.0; R=0.0
    for j in range(0,m):
        nj=layers[j][0]
        Dj=layers[j][1]
        (r,t,s)=rts(n0, nj, Dj)
        Vlast=V
        (U,V)=(r*V + s*U, V - r*U)
        T=T*t*Vlast/V
        R=U/V
    return (R,T)

def logTm(m, n0, layers):
    "log of transmission coefficient for m layers"
    U=0.0; V=1.0; logT=0.0
    for j in range(0,m):
        nj=layers[j][0]
        Dj=layers[j][1]
        (r,t,s)=rts(n0, nj, Dj)
        (U,V)=(r*V + s*U, V - r*U)
        logT=logT + logT1(n0, nj, Dj)

    return logT - numpy.log(V)

def RTmdash(m, n0, layers, q):
    "Derivative of R and T for m layers wrt n_q"
    U=0.0; V=1.0; T=1.0; R=0.0
    Udash=0.0; Vdash=0.0; Tdash=0.0
    for j in range(0,m):
        nj=layers[j][0]
        Dj=layers[j][1]
        (r,t,s)=rts(n0, nj, Dj)
        Vlast=V; Vdashlast=Vdash
        if j==q: (rdash,tdash,sdash)=rtsdash(n0, nj, Dj)
        else:    (rdash,tdash,sdash)=(0.0,0.0,0.0)
        (U, V, Udash, Vdash)=(r*V + s*U, V - r*U, \
            rdash*V + r*Vdash + sdash*U + s*Udash, Vdash - rdash*U - r*Udash)
        Tlast=T; Tdashlast=Tdash
        T=T*t*Vlast/V
        Tdash=T*(Tdashlast/Tlast + tdash/t + Vdashlast/Vlast - Vdash/V)
    return Tdash

def rts(n0, nj, Dj):
    "returns tuple of r, t and t**2 - r**2 for layer j"
    c=numpy.cos(nj*Dj)
    s=numpy.sin(nj*Dj)
    if nj!=0:
        d=c + (0.5J)*(nj/n0 + n0/nj)*s
        r=(0.5J)*s*(n0/nj - nj/n0)/d
        t=1.0/d
    else:
        d=r=t=0
    return (r, t*numpy.exp((1J)*n0*Dj), t*t-r*r)

def rtsdash(n0, nj, Dj):
    "Derivative of r and t for layer j"
    c=numpy.cos(nj*Dj)
    s=numpy.sin(nj*Dj)
    (r,t,junk)=rts(n0, nj, Dj)
    d=c + (0.5J)*(nj/n0 + n0/nj)*s
    ddash= (-Dj + (0.5J)*(1.0/n0 - n0/nj**2))*s + (0.5J)*Dj*(nj/n0 + n0/nj)*c
    rdash= r*((nj**2+n0**2)/(nj*(nj**2-n0**2)) + Dj*c/s - ddash/d)
    tdash= -t*ddash/d
    sdash= 2.0*(tdash*t - rdash*r)
    return (rdash, tdash, sdash)


############################################################################
# Functions for creating simulated data
############################################################################

def simulate_parallel(reference_file, layers, dk, noise_level):
    x=ops.readdatafile(reference_file, 1.0)
    refindex=ops.readdatafile2("i500infR.txt", 1.0)
    refextin=ops.readdatafile2("e500infR.txt", 1.0)
    N=len(x)
    m=len(layers)
    indices=numpy.array([l[0] for l in layers])
    thicknesses=numpy.array([l[1] for l in layers])
    T=numpy.zeros(N, complex)
    M=N/2
    for k in range(0,M):
        kD=k*dk*thicknesses
        'use a constant complex refractive index'
        nkD=numpy.array([(indices[l],kD[l]) for l in range(0,m)])
        'import complex refractive index from refindex and refextin'
        #nkD=numpy.array([(refindex[k]+(1.0J)*refextin[k],kD[l]) for l in range(0,m)])
        T[k]=RTm(m, 1.0, nkD)[1]
        print T[k]
    T[0]=0.0+0.0J
    for k in range(1,M): 
        T[N-k]=numpy.conjugate(T[k])
    X=fft(x)/N
    Y=T*X
    g=numpy.zeros(N, float)
    for t in range(0,N): g[t]=random.gauss(0.00, noise_level)
    return N*ifft(Y).real + g

def simulate_noR(n0, nkD):
    n1=nkD[0][0]
    kD=nkD[0][1]
    Tab=2.0*n1/(n1+n0)
    Tba=2.0*n0/(n1+n0)
    Rab=(n1-n0)/(n1+n0)
    Pair=numpy.exp(-(1.0J)*n0*kD)
    Psam=numpy.exp(-(1.0J)*n1*kD)
    if n1!=0:
        a=Tab*Psam*Tba/Pair
    else:
        a=0
    return a

def simulate_infR(n0, nkD):
    n1=nkD[0][0]
    kD=nkD[0][1]   
    Tab=2.0*n1/(n1+n0)
    Tba=2.0*n0/(n1+n0)
    Rab=(n1-n0)/(n1+n0)
    Pair=numpy.exp(-(1.0J)*n0*kD)
    Psam=numpy.exp(-(1.0J)*n1*kD)
    if n1!=0:
        a=Tab*Psam*Tba
        b=Pair*(1.0-Rab*Rab*Psam*Psam)
    else:
        a=0
        b=1
    return a/b

def simulate_oneR(n0, nkD):
    n1=nkD[0][0]
    kD=nkD[0][1]
    Tab=2.0*n1/(n1+n0)
    Tba=2.0*n0/(n1+n0)
    Rab=(n1-n0)/(n1+n0)
    Pair=numpy.exp(-(1.0J)*n0*kD)
    Psam=numpy.exp(-(1.0J)*n1*kD)
    if n1!=0:
        a=Tab*Psam*Tba/Pair*(1.0+Rab*Rab*Psam*Psam)
    else:
        a=0
    return a

def simulate_parallel2(reference_file, layers, dk, noise_level, paddingfactor, N, M, numRefl_simulated_plane):
    xx=ops.readdatafile(reference_file, 1.0, paddingfactor)
    x=numpy.zeros(N,float)
    x[0:M]=xx[0:M]
    #refindex=ops.readdatafile2("i500infR.txt", 1.0)
    #refextin=ops.readdatafile2("e500infR.txt", 1.0)
   
    #custom ext coefficient#
    #def exttoabs(ni,dk,k):
    #    return 2*numpy.pi*dk*ni*k

    #def abstoext(abs,dk,k):
    #    return abs/2*numpy.pi*dk*k

    #exttoabs2=numpy.zeros(2048)
    #for i in numpy.arange(0,2048):
    #    exttoabs2[i]=exttoabs(0.005,dk,i)
    
    #abstoext2=numpy.zeros(2048)
    #for i in numpy.arange(0,2048):
    #    abstoext2[i]=abstoext(0.00005,dk,i)
    
    
    #refindex=numpy.zeros(2048)
    #for i in numpy.arange(0,2048):
    #    refindex[i]=2.000
        
    #refextin=abstoext2
    #refextin=abstoext2
    
    #pylab.figure()
    #pylab.plot(fscale[P:Q], exttoabs2[P:Q])
    #pylab.title('abs')
    #pylab.legend()
    #pylab.grid(True)  

    #pylab.figure()
    #pylab.plot(fscale[P:Q], abstoext2[P:Q])
    #pylab.title('ext')
    #pylab.legend()
    #pylab.grid(True)  
    
    
    m=len(layers)
    indices=numpy.array([l[0] for l in layers])
    thicknesses=numpy.array([l[1] for l in layers])
    T=numpy.zeros(N, complex)
    Nhalf=N/2
    for k in range(0,Nhalf):
        kD=k*dk*thicknesses
        'use a constant complex refractive index'
        nkD=numpy.array([(indices[l],kD[l]) for l in range(0,m)])
        'import complex refractive index from refindex and refextin'
        #nkD=numpy.array([(refindex[k]+(1.0J)*refextin[k],kD[l]) for l in range(0,m)])
        
        if numRefl_simulated_plane == 'none':
            T[k]=simulate_noR(1.0,nkD)

        if  numRefl_simulated_plane == 'inf':
            T[k]=simulate_noR(1.0,nkD)
   
        if  numRefl_simulated_plane == 'one':
            T[k]=simulate_noR(1.0,nkD)     

    T[0]=0.0+0.0J
    for k in range(1,Nhalf): 
        T[N-k]=numpy.conjugate(T[k])
    X=fft(x)/N
    Y=T*X
    g=numpy.zeros(N, float)
    for t in range(0,N): g[t]=random.gauss(0.00, noise_level)
    return N*ifft(Y).real + g

def simulate_converging_freq(reference_file, n0, n1, thickness, theta_cutoff, dk, noise_level, paddingfactor, N, M, polarisation_configuration_simulated, numRefl_simulated):
    xx=ops.readdatafile(reference_file, 1.0, paddingfactor)
    x=numpy.zeros(N,float)  
    x[0:M]=xx[0:M]  
    truncateval=50   
    x=ops.truncate(x, truncateval, 30)    
    T=converging_freq(n0, n1, thickness, theta_cutoff, dk, N, polarisation_configuration_simulated, numRefl_simulated)     
    X=fft(x)/N
    Y=T*X    
    g=numpy.zeros(N, float)
    for t in range(0,N): g[t]=random.gauss(0.0, noise_level)
    return N*ifft(Y).real + g

def converging_freq(n0, n1, thickness, angular_spread, dk, N, polarisation_configuration_simulated, numRefl_simulated):
    M=N/2
    TF=numpy.zeros(N,complex)

    i=0
    for k in range(0,M):        
        #set up arrays of thetas to be used in integration
        # theta_max is maximum value of theta used in numerical integration
        theta_cutoff=angular_spread[k]
        theta_max=2.0*theta_cutoff
        theta=0.5*theta_max*(1.0 + xi)      
        weight=0.5*wi*theta_max*angular_cutoff_simulated(theta, theta_cutoff)
        #1.0*numpy.exp(-1.0*theta**2.0/theta0**2.0)
        # compute reflection coeffs for these angles
        # these are independent of freq
        theta_i=theta
        theta_t=numpy.arcsin(n0*numpy.sin(theta_i)/n1[k])
        (R1_para, R1_perp, T1_para, T1_perp)=fresnel_oblique(n0, n1[k], theta_i, theta_t)
        (R2_para, R2_perp, T2_para, T2_perp)=fresnel_oblique(n1[k], n0, theta_t, theta_i)
        #
        T_para = T1_para * T2_para
        R_para = R2_para * R2_para
        T_perp = T1_perp * T2_perp
        R_perp = R2_perp * R2_perp
        s=numpy.sin(theta)
        c=numpy.cos(theta)
        kD=k*dk*thickness
# P1, P2 are arrays over theta for this freq
        P1=numpy.exp(-1.0J*kD*(n1[k]*numpy.cos(theta_t) - n0*numpy.cos(theta_i)))
        P2=numpy.exp(-1.0J*kD*n1[k]*numpy.cos(theta_t)) 
               
        if numRefl_simulated == 'none':

            F_para=T_para*P1
            F_perp=T_perp*P1
 
        if  numRefl_simulated == 'inf':

            F_para=T_para*P1/(1.0 - P2*P2*R_para)
            F_perp=T_perp*P1/(1.0 - P2*P2*R_perp)      

        if  numRefl_simulated == 'one':

            F_para=T_para*P1*(1.0 + P2*P2*R_para)
            F_perp=T_perp*P1*(1.0 + P2*P2*R_perp)
        
        if  polarisation_configuration_simulated == 'PM':
            "polarised magnetic field configuration"
            if theta_cutoff == 0:
                TF[i] = 0
            else:
                TF[i]=sum((c*c*F_para + F_perp)*s*weight)/sum((1.0 + c*c)*s*weight)
                       
        if  polarisation_configuration_simulated == 'PE':
            "polarised electric field configuration"
            if theta_cutoff == 0:  
                TF[i] = 0
            else:
                TF[i]=sum((F_para + F_perp)*c*s*weight)/sum(2.0*c*s*weight)      
        i=i+1  
                         
    for i in range(1,M):
        TF[N-i]=numpy.conjugate(TF[i])
  
    
    return TF

def simulate_converging(reference_file, n0, n1, thickness, theta_cutoff, dk, noise_level, paddingfactor, N, M, polarisation_configuration_simulated, numRefl_simulated):
    xx=ops.readdatafile(reference_file, 1.0, paddingfactor)
    x=numpy.zeros(N,float)
    x[0:M]=xx[0:M]    
    T=converging(n0, n1, thickness, theta_cutoff, dk, N, polarisation_configuration_simulated, numRefl_simulated)
    X=fft(x)/N
    Y=T*X
    g=numpy.zeros(N, float)
    for t in range(0,N): g[t]=random.gauss(0.0, noise_level)
    return N*ifft(Y).real + g

def simulate_converging_conical(reference_file, n0, n1, thickness, theta_c, dk, noise_level, paddingfactor, N, M):
    xx=ops.readdatafile(reference_file, 1.0, paddingfactor)
    x=numpy.zeros(N,float)
    x[0:M]=xx[0:M]
    T=converging_conical(n0, n1, thickness, theta_c, dk, N)
    X=fft(x)/N
    Y=T*X
    g=numpy.zeros(N, float)
    for t in range(0,N): g[t]=random.gauss(0.0, noise_level)
    return N*ifft(Y).real + g


############################################################################
# Functions for using a converging beam
############################################################################
#
# Points and weights for gaussian integration
xi16=[0.09501250983764, 0.28160355077926, 0.45801677765723, 0.61787624440264, \
      0.75540440835500, 0.86563120238783, 0.94457502307323, 0.98940093499165]
wi16=[0.18945061045507, 0.18260341504492, 0.16915651939500, 0.14959598881658, \
      0.12462897125553, 0.09515851168249, 0.06225352393865, 0.02715245941175]

xi32=[0.04830766568774, 0.14447196158280, 0.23928736225214, 0.33186860228213, \
      0.42135127613064, 0.50689990893223, 0.58771575724076, 0.66304426693022, \
      0.73218211874029, 0.79448379596794, 0.84936761373257, 0.89632115576605, \
      0.93490607593774, 0.96476225558751, 0.98561151154527, 0.99726386184948]
wi32=[0.09654008851473, 0.09563872007927, 0.09384439908080, 0.09117387869576, \
      0.08765209300440, 0.08331192422695, 0.07819389578707, 0.07234579410885, \
      0.06582222277636, 0.05868409347854, 0.05099805926238, 0.04283589802223, \
      0.03427386291302, 0.02539206530926, 0.01627439473091, 0.00701861000947]

# use 16 point gauss formula, should be accurate enough
#xi=numpy.array(xi16+list(-numpy.array(xi16)))
#wi=numpy.array(wi16+wi16)

# use 32 point gauss formula, slower but more accurate
iw32=wi32[:]
iw32.reverse()
ix32=xi32[:]
ix32.reverse()

ix=numpy.array(ix32+list(-numpy.array(xi32)))
xi=numpy.array(xi32+list(-numpy.array(xi32)))
#xi=numpy.array(xi32+list(numpy.array(xi32)-1))
wi=numpy.array(wi32+wi32)
iw=numpy.array(wi32+list(-numpy.array(wi32)))

def converging_conical(n0, n1, thickness, theta_c, dk, N):
    M=N/2
#set up arrays of thetas to be used in integration
# theta_max is maximum value of theta used in numerical integration
    theta_w=0.01
    theta_max=2.0*theta_w
    theta=theta_max*xi+theta_c
    angular_cutoff=numpy.exp(-2.0*(theta-theta_c)**2.0/theta_w**2.0)   
    "donut cutoff with no f dependence (e-2) theta0 is the half angle"
    weight=wi*theta_max*angular_cutoff
# compute reflection coeffs for these angles
# these are independent of freq
    theta_i=theta
    theta_t=numpy.arcsin(n0*numpy.sin(theta_i)/n1)
    (R1_para, R1_perp, T1_para, T1_perp)=fresnel_oblique(n0, n1, theta_i, theta_t)
    (R2_para, R2_perp, T2_para, T2_perp)=fresnel_oblique(n1, n0, theta_t, theta_i)
#
    T_para = T1_para * T2_para
    R_para = R2_para * R2_para
    T_perp = T1_perp * T2_perp
    R_perp = R2_perp * R2_perp
    s=numpy.sin(theta)
    c=numpy.cos(theta)
#
    TF=numpy.zeros(N,complex)
    i=0
    for k in range(0,M+1):
        kD=k*dk*thickness
# P1, P2 are arrays over theta for this freq
        P2=numpy.exp(-1.0J*kD*n1/numpy.cos(theta_t))
#        P1=numpy.exp(-1.0J*kD*n0/numpy.cos(theta0))
        P3=numpy.exp(-1.0J*kD*(n1*numpy.cos(theta_t) - n0*numpy.cos(theta_i)))
#        F_para=T_para*(P2/P1)/(1.0 - P2*P2*R_para)
#        F_perp=T_perp*(P2/P1)/(1.0 - P2*P2*R_perp)
        "inf reflections"
        #F_para=T_para*P3/(1.0 - P2*P2*R_para)
        #F_perp=T_perp*P3/(1.0 - P2*P2*R_perp)
        "zero reflections"
        F_para=T_para*P3
        F_perp=T_perp*P3
        TF[i]=sum((c*c*F_para + F_perp)*s*weight)
        i=i+1
    for i in range(1,M): TF[N-i]=numpy.conjugate(TF[i])
    norm=sum((1.0 + c*c)*s*weight)
    return TF/norm

def converging(n0, n1, thickness, angular_spread, dk, N, polarisation_configuration_simulated, numRefl_simulated):
    M=N/2
    TF=numpy.zeros(N,complex)
    i=0
    for k in range(0,M):        
        #set up arrays of thetas to be used in integration
        # theta_max is maximum value of theta used in numerical integration
        theta_cutoff=angular_spread[k]
        theta_max=2.0*theta_cutoff
        theta=0.5*theta_max*(1.0 + xi)
        weight=0.5*wi*theta_max*angular_cutoff_simulated(theta, theta_cutoff)
        # compute reflection coeffs for these angles
        # these are independent of freq
        theta_i=theta
        theta_t=numpy.arcsin(n0*numpy.sin(theta_i)/n1)
        (R1_para, R1_perp, T1_para, T1_perp)=fresnel_oblique(n0, n1, theta_i, theta_t)
        (R2_para, R2_perp, T2_para, T2_perp)=fresnel_oblique(n1, n0, theta_t, theta_i)
        #
        T_para = T1_para * T2_para
        R_para = R2_para * R2_para
        T_perp = T1_perp * T2_perp
        R_perp = R2_perp * R2_perp
        s=numpy.sin(theta)
        c=numpy.cos(theta)
        #
        kD=k*dk*thickness
# P1, P2 are arrays over theta for this freq
        P1=numpy.exp(-1.0J*kD*(n1*numpy.cos(theta_t) - n0*numpy.cos(theta_i)))
        P2=numpy.exp(-1.0J*kD*n1*numpy.cos(theta_t))   
        
        if numRefl_simulated == 'none':

            F_para=T_para*P1
            F_perp=T_perp*P1

        if  numRefl_simulated == 'inf':

            F_para=T_para*P1/(1.0 - P2*P2*R_para)
            F_perp=T_perp*P1/(1.0 - P2*P2*R_perp)

        if  numRefl_simulated == 'one':

            F_para=T_para*P1*(1.0 + P2*P2*R_para)
            F_perp=T_perp*P1*(1.0 + P2*P2*R_perp)
        
        if  polarisation_configuration_simulated == 'PM':
            "polarised magnetic field configuration"
            TF[i]=sum((c*c*F_para + F_perp)*s*weight)

        
        if  polarisation_configuration_simulated == 'PE':
            "polarised electric field configuration"
            TF[i]=sum((F_para + F_perp)*c*s*weight)
           
        i=i+1  
                         
    for i in range(1,M): TF[N-i]=numpy.conjugate(TF[i])

    if  polarisation_configuration_simulated == 'PM':
        "polarised magnetic field configuration"
        norm=sum((1.0 + c*c)*s*weight)    
    
    if  polarisation_configuration_simulated == 'PE':
        "polarised electric field configuration"
        norm=sum(2.0*c*s*weight)  
    
    return TF/norm
        
    "for the simulated data"
def angular_cutoff_simulated(theta, theta0):
    "model for cutting off angular spectrum at theta0"
    #return 1.0/(1.0 + (theta/theta0)**6)
    
    "gaussian cutoff with no f dependence (e-1) theta0 is the half angle"
    return 1.0*numpy.exp(-1.0*theta**2.0/theta0**2.0)
    "gaussian cutoff with no f dependence (e-2) theta0 is the half angle"
    #return 1.0*numpy.exp(-2.0*theta**2.0/theta0**2.0)
    "square cutoff theta_max determines 2*width"
    #return theta/theta
    "square cutoff where theta0 determines the halfwidth"
    #return theta0  



def fresnel_oblique(ni, nt, theta_i, theta_t):
    "Fresnel reflection and transmission coefficients for oblique incidence"
    ci=numpy.cos(theta_i)
    ct=numpy.cos(theta_t)
    d_perp=ni*ci + nt*ct
    d_para=ni*ct + nt*ci
    r_perp=(ni*ci-nt*ct)/d_perp
    r_para=(nt*ci-ni*ct)/d_para
    t_perp=2.0*ni*ci/d_perp
    t_para=2.0*ni*ci/d_para
    return (r_para, r_perp, t_para, t_perp)

def angular_cutoff(theta, theta0, freqk, kmin):
    "model for cutting off angular spectrum at theta0"
    "original"
    #return 1.0/(1.0 + (theta*(1 + float(freqk)/kmin)/theta0)**6)
    "no f dependence"
    #return 1.0/(1.0 + (theta/theta0)**6)
    "new function"
    #return 1.0/(1.0 + ((theta/theta0)**6)*(1.0+(float(freqk)/kmin)**6))
    "gaussian cutoff with no f dependence (FWHM)"
    #rho=theta0/(2.0*numpy.sqrt(2.0*numpy.log(2.0)))
    #return 1.0*numpy.exp(-theta**2.0/(2.0*rho**2.0))
    "square cutoff theta_max determines 2*width"
    #return theta/theta
    "square cutoff where theta0 determines the halfwidth"
    #return theta0
    "gaussian cutoff with f dependence"
    #rho=theta0/(2.0*numpy.sqrt(2.0*numpy.log(2.0)))*(1.0+(freqk/float(kmin))**6.0)
    #return 1.0*numpy.exp(-theta**2.0/(2.0*rho**2.0))
    "gaussian cutoff with width dependence"
    #rho=1.0/(2.0*numpy.sqrt(2.0*numpy.log(2.0))*dk*freqk*width2[freqk])
    #return 1.0*numpy.exp(-theta**2.0/(2.0*rho**2.0))
    "gaussian cutoff with no f dependence (e-2) theta0 is the half angle"
    #return 1.0*numpy.exp(-2.0*theta**2.0/theta0**2.0)
    "donut cutoff with no f dependence (e-2) theta0 is the half angle"
    #return 1.0*numpy.exp(-2*(theta-0.2)**2.0/theta0**2.0)
    
    "gaussian cutoff with no f dependence (e-2) theta0 is the full angle"
    return 1.0*numpy.exp(-1.0*theta**2.0/theta0**2.0)
    
    
def Tconv(n0, n1, kD, freqk, polarisation_configuration, numRefl_conv):
    "fitting function for converging beam"
    set_up_converging_fit(theta_cutoff, freqk, polarisation_configuration)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    sinplus = sintheta*costheta_t + costheta*sintheta_t
    cosplus = costheta*costheta_t - sintheta*sintheta_t
    sinminus= sintheta*costheta_t - costheta*sintheta_t
    cosminus= costheta*costheta_t + sintheta*sintheta_t
    r_perpsq=(sinminus/sinplus)**2
    r_parasq=r_perpsq*(cosplus/cosminus)**2
    
    P1=numpy.exp(-1.0J*kD*(n1*costheta_t - n0*costheta))
    
    P2=numpy.exp(-1.0J*kD*n1*costheta_t)
    t_perp12=4.0*costheta*sintheta_t*costheta_t*sintheta/sinplus**2
    t_para12=t_perp12/cosminus**2
    
    if numRefl_conv == 'none':
        F_para=t_para12*P1
        F_perp=t_perp12*P1

    if  numRefl_conv == 'inf': 
        F_para=t_para12*P1/(1.0 - P2*P2*r_parasq)
        F_perp=t_perp12*P1/(1.0 - P2*P2*r_perpsq)

    if  numRefl_conv == 'one':
        F_para=t_para12*P1*(1.0 + P2*P2*r_parasq)
        F_perp=t_perp12*P1*(1.0 + P2*P2*r_perpsq)    
    
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        return sum((costheta*costheta*F_para + F_perp)*sintheta*weight)/fit_norm

    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        return sum((F_para + F_perp)*costheta*sintheta*weight)/fit_norm   


def set_up_converging_fit(theta_cutoff, freqk, polarisation_configuration):
    "creates arrays theta and weight for fitting function Tconv"
    global theta, weight, sintheta, costheta, fit_norm
    theta_max=2.0*theta_cutoff
    theta=0.5*theta_max*(1.0 + xi)
    #weight=0.5*wi*theta_max*theta_cutoff
    weight=0.5*wi*theta_max*angular_cutoff(theta, theta_cutoff, freqk, kmin)
    profile=angular_cutoff(theta, theta_cutoff, freqk, kmin)
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        fit_norm=sum((1.0 + costheta*costheta)*sintheta*weight)

    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        fit_norm=sum(2.0*costheta*sintheta*weight)
        
def set_theta_cutoff(theta, freqk_min, shiftguassian):
    "used to get parameters into module 'models'"
    global theta_cutoff, kmin, angularcentre
    theta_cutoff=theta
    kmin=freqk_min
    angularcentre=shiftguassian
    
#theta_fdcutoff=numpy.zeros(N)
#theta_fdcutoff[:70]=0.28
#theta_fdcutoff[70:]=0.24

#deltat=1.09e-14
#deltaf=1.0/(N*deltat)

#tscale=numpy.zeros(N, complex)
#fscale=numpy.zeros(N, complex)
#for i in range(0, N):
#    tscale[i]=i*deltat*1e12
#    fscale[i]=i*deltaf*1e-12
#pylab.figure()
#pylab.plot(fscale[0:512], theta_fdcutoff[0:512])
#pylab.title('fdcutoff')
#pylab.legend()
#pylab.grid(True)
#pylab.show()

def angular_cutoffdonut(theta, theta0, freqk, angularcentre):
    return 1.0*numpy.exp(-2*(theta-angularcentre)**2.0/theta0**2.0)    

def set_up_converging_fitdonut(theta_cutoff, freqk):
    "creates arrays theta and weight for fitting function Tconv"
    global theta, weight, sintheta, costheta, fit_norm
    theta_max=2.0*theta_cutoff
    theta=theta_max*xi+angularcentre
    "donut cutoff with no f dependence (e-2) theta0 is the half angle"
    weight=wi*theta_max*angular_cutoffdonut(theta, theta_cutoff, freqk, angularcentre)
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    fit_norm=sum((1.0 + costheta*costheta)*sintheta*weight)

def Tconvdonut(n0, n1, kD, freqk):
    "fitting function for converging beam"
    set_up_converging_fitdonut(theta_cutoff, freqk)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    sinplus = sintheta*costheta_t + costheta*sintheta_t
    cosplus = costheta*costheta_t - sintheta*sintheta_t
    sinminus= sintheta*costheta_t - costheta*sintheta_t
    cosminus= costheta*costheta_t + sintheta*sintheta_t
    r_perpsq=(sinminus/sinplus)**2
    r_parasq=r_perpsq*(cosplus/cosminus)**2
    P3=numpy.exp(-1.0J*kD*(n1*costheta_t - n0*costheta))
    P2=numpy.exp(-1.0J*kD*n1/costheta_t)
    #P1=exp(-1.0J*kD*n0/costheta)
    t_perp12=4.0*costheta*sintheta_t*costheta_t*sintheta/sinplus**2
    t_para12=t_perp12/cosminus**2
    "inf Reflections"
    #F_para=t_para12*P3/(1.0 - P2*P2*r_parasq)
    #F_perp=t_perp12*P3/(1.0 - P2*P2*r_perpsq)
    "zero Reflection"
    F_para=t_para12*P3
    F_perp=t_perp12*P3
    "1 Reflection"
    #F_para=t_para12*P3*(1.0+P2*P2*r_parasq)
    #F_perp=t_perp12*P3*(1.0+P2*P2*r_perpsq)
    return sum((costheta*costheta*F_para + F_perp)*sintheta*weight)/fit_norm




    
def angular_cutoff_theta_fdcutoff(theta, theta_0, freqk, kmin):
    "gaussian cutoff with no f dependence (e-2) theta0 is the half angle"
    #return 1.0*numpy.exp(-2.0*theta**2.0/theta_0**2.0)
    "gaussian cutoff with no f dependence (e-2) theta0 is the full angle"
    return 1.0*numpy.exp(-1.0*theta**2.0/theta_0**2.0)
    
def set_up_converging_fit_theta_fdcutoff(theta_0, freqk, polarisation_configuration):
    "creates arrays theta and weight for fitting function Tconv"
    global theta, weight, sintheta, costheta, fit_norm
    theta_max=2.0*theta_0
    theta=0.5*theta_max*(1.0 + xi)
    #weight=0.5*wi*theta_max*theta_cutoff

    weight=0.5*wi*theta_max*angular_cutoff_theta(theta, theta_0, freqk, kmin)
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        fit_norm=sum((1.0 + costheta*costheta)*sintheta*weight)     
        
    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        fit_norm=sum(2.0*costheta*sintheta*weight)        
                
def Tconvfd(n0, n1, kD, freqk, theta_fdcutoff, polarisation_configuration):
    "fitting function for converging beam"     
    set_up_converging_fit_theta_fdcutoff(theta_fdcutoff, freqk, polarisation_configuration)
       
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    sinplus = sintheta*costheta_t + costheta*sintheta_t
    cosplus = costheta*costheta_t - sintheta*sintheta_t
    sinminus= sintheta*costheta_t - costheta*sintheta_t
    cosminus= costheta*costheta_t + sintheta*sintheta_t
    r_perpsq=(sinminus/sinplus)**2
    r_parasq=r_perpsq*(cosplus/cosminus)**2
    P3=numpy.exp(-1.0J*kD*(n1*costheta_t - n0*costheta))
    P2=numpy.exp(-1.0J*kD*n1*costheta_t)
    #P1=exp(-1.0J*kD*n0/costheta)
    t_perp12=4.0*costheta*sintheta_t*costheta_t*sintheta/sinplus**2
    t_para12=t_perp12/cosminus**2
    "inf Reflections"
    #F_para=t_para12*P3/(1.0 - P2*P2*r_parasq)
    #F_perp=t_perp12*P3/(1.0 - P2*P2*r_perpsq)
    "zero Reflection"
    F_para=t_para12*P3
    F_perp=t_perp12*P3
    "1 Reflection (this is probably wrong needs checking"
    #F_para=t_para12*P3*(1.0+P2*P2*r_parasq)
    #F_perp=t_perp12*P3*(1.0+P2*P2*r_perpsq)
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        return sum((costheta*costheta*F_para + F_perp)*sintheta*weight)/fit_norm

    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        return sum((F_para + F_perp)*costheta*sintheta*weight)/fit_norm
    
    
"Angle dependent converging beam transfer functions"


def angular_cutoff_thetadonut(theta, theta_w, freqk, theta_c):
    return 1.0*numpy.exp(-2.0*(theta-theta_c)**2.0/theta_w**2.0)    

def set_up_converging_fit_thetadonut(theta_c, freqk, polarisation_configuration):
    "creates arrays theta and weight for fitting function Tconv"
    global theta, weight, sintheta, costheta, fit_norm
    theta_w=0.01
    theta_max=2.0*theta_w
    theta=theta_max*xi+theta_c
    "donut cutoff with no f dependence (e-2) theta0 is the half angle"
    weight=wi*theta_max*angular_cutoff_thetadonut(theta, theta_w, freqk, theta_c)
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        fit_norm=sum((1.0 + costheta*costheta)*sintheta*weight)

    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        fit_norm=sum(2.0*costheta*sintheta*weight)
    


def angular_cutoff_theta(theta, theta_0, freqk, kmin):
    "gaussian cutoff with no f dependence (e-2) theta0 is the half angle"
    #return 1.0*numpy.exp(-2.0*theta**2.0/theta_0**2.0)
    "gaussian cutoff with no f dependence (e-2) theta0 is the full angle"
    return 1.0*numpy.exp(-1.0*theta**2.0/theta_0**2.0)
    "gaussian cutoff with no f dependence (e-2) theta0 is the full angle"
    #return 2.0*numpy.exp(-1.0*theta**2.0/theta_0**2.0)/theta_0**2
    
def set_up_converging_fit_theta(theta_0, freqk, polarisation_configuration):
    "creates arrays theta and weight for fitting function Tconv"
    global theta, weight, sintheta, costheta, fit_norm
    theta_max=2.0*theta_0
    theta=0.5*theta_max*(1.0 + xi)
    #weight=0.5*wi*theta_max*theta_cutoff

    weight=0.5*wi*theta_max*angular_cutoff_theta(theta, theta_0, freqk, kmin)
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        fit_norm=sum((1.0 + costheta*costheta)*sintheta*weight)     
        
    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        fit_norm=sum(2.0*costheta*sintheta*weight)        
        
def Tconvtheta(n0, n1, kD, freqk, theta_0, polarisation_configuration):
    "fitting function for converging beam"

    set_up_converging_fit_theta(theta_0, freqk, polarisation_configuration)
    #set_up_converging_fit_thetadonut(theta_0, freqk,  polarisation_configuration)
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    sinplus = sintheta*costheta_t + costheta*sintheta_t
    cosplus = costheta*costheta_t - sintheta*sintheta_t
    sinminus= sintheta*costheta_t - costheta*sintheta_t
    cosminus= costheta*costheta_t + sintheta*sintheta_t
    r_perpsq=(sinminus/sinplus)**2
    r_parasq=r_perpsq*(cosplus/cosminus)**2
    P3=numpy.exp(-1.0J*kD*(n1*costheta_t - n0*costheta))
    P2=numpy.exp(-1.0J*kD*n1*costheta_t)
    #P1=exp(-1.0J*kD*n0/costheta)
    t_perp12=4.0*costheta*sintheta_t*costheta_t*sintheta/sinplus**2
    t_para12=t_perp12/cosminus**2
    "inf Reflections"
    #F_para=t_para12*P3/(1.0 - P2*P2*r_parasq)
    #F_perp=t_perp12*P3/(1.0 - P2*P2*r_perpsq)
    "zero Reflection"
    F_para=t_para12*P3
    F_perp=t_perp12*P3
    "1 Reflection"
    #F_para=t_para12*P3*(1.0+P2*P2*r_parasq)
    #F_perp=t_perp12*P3*(1.0+P2*P2*r_perpsq)
  
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        return sum((costheta*costheta*F_para + F_perp)*sintheta*weight)/fit_norm

    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        return sum((F_para + F_perp)*costheta*sintheta*weight)/fit_norm




def set_up_converging_fit_theta_cone(theta_0, freqk, polarisation_configuration):
    "creates arrays theta and weight for fitting function Tconv"
    global theta, weight, sintheta, costheta, fit_norm
    theta=theta_0
    weight=1.0
    sintheta=numpy.sin(theta)
    costheta=numpy.cos(theta)
    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        fit_norm=((1.0 + costheta*costheta)*sintheta*weight)

    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        fit_norm=(2.0*costheta*sintheta*weight)        

def Tconvthetacone(n0, n1, kD, freqk, theta_0, polarisation_configuration):
    
    "fitting function for converging beam"
    
    set_up_converging_fit_theta_cone(theta_0, freqk, polarisation_configuration)
    
    theta_t=numpy.arcsin(sintheta/n1)
    costheta_t=numpy.cos(theta_t)
    sintheta_t=numpy.sin(theta_t)
    sinplus = sintheta*costheta_t + costheta*sintheta_t
    cosplus = costheta*costheta_t - sintheta*sintheta_t
    sinminus= sintheta*costheta_t - costheta*sintheta_t
    cosminus= costheta*costheta_t + sintheta*sintheta_t
    r_perpsq=(sinminus/sinplus)**2
    r_parasq=r_perpsq*(cosplus/cosminus)**2
    P3=numpy.exp(-1.0J*kD*(n1*costheta_t - n0*costheta))
    P2=numpy.exp(-1.0J*kD*n1*costheta_t)
    #P1=exp(-1.0J*kD*n0/costheta)
    t_perp12=4.0*costheta*sintheta_t*costheta_t*sintheta/sinplus**2
    t_para12=t_perp12/cosminus**2
    "inf Reflections"
    #F_para=t_para12*P3/(1.0 - P2*P2*r_parasq)
    #F_perp=t_perp12*P3/(1.0 - P2*P2*r_perpsq)
    "zero Reflection"
    F_para=t_para12*P3
    F_perp=t_perp12*P3
    "1 Reflection"
    #F_para=t_para12*P3*(1.0+P2*P2*r_parasq)
    #F_perp=t_perp12*P3*(1.0+P2*P2*r_perpsq)

    if  polarisation_configuration == 'PM':
        "polarised magnetic field configuration"
        return ((costheta*costheta*F_para + F_perp)*sintheta*weight)/fit_norm

    if  polarisation_configuration == 'PE':
        "polarised electric field configuration"
        return ((F_para + F_perp)*costheta*sintheta*weight)/fit_norm



#def conv_ri(kstart, thickscale, model, deriv):
#    global T0, kD
# start at good freq and step down
#    for k in range(kstart, 0, -1):
#        T0=Y[k]/X[k]
#        kD=thickscale*k*kl
#        z0=n[k]
#        (ok,z0)=newton(z0, model, deriv)
#        if not ok: break
#        n[k]=z0
#    k_lo =k
# start at kstart and step up
#    for k in range(kstart, M):
#        T0=Y[k]/X[k]
#        kD=thickscale*k*kl
#        z0=n[k]
#        (ok,z0)=newton(z0, model, deriv)
#        if not ok: break
#        n[k]=z0
#    k_hi=k
#    return (k_lo, k_hi, n)
