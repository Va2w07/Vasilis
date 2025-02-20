from numpy import *
from numpy.fft import fft, ifft
#from cmath import *      # note this prevents Ufuncs existing
import ops
from physics import graph, clear_graph, pause

global T0, kD
T0=0.0; kD=0.0

def set_thickscale(s):
    global thickscale
    thickscale=s

def set_nbread(n):
    global nbread
    nbread=n

def set_nfill(n):
    global nfill
    nfill=n

##############################################################################
# Models and derivatives fitted to transfer function
##############################################################################

# Single slab of refractive index n in free space
# fitting raw transfer function, used only for refinement
def rawslab1(n):
    return T1(1.0, n, thickscale*ops.kD[0]) - ops.T0

def rawslab1dash(n):
    return T1dash(1.0, n, thickscale*ops.kD[0])

# fitting log of transfer function, normal single layer function
def slab1(n):
    return logT1(1.0, n, thickscale*ops.kD[0]) - ops.T0

def slab1dash(n):
    return logT1dash(1.0, n, thickscale*ops.kD[0])

# Three layer sandwich, fitting middle layer
def GBsand(n):
    return logT3(1.0, nbread, n, nbread, \
                 ops.kD[0], thickscale*ops.kD[1], ops.kD[2]) - ops.T0

def GBsanddash(n):
    return logT3dash(1.0, nbread, n, nbread, \
                     ops.kD[0], thickscale*ops.kD[1], ops.kD[2])

# Filled roll, fitting outer layers
def filled_roll(n):
    return logT3(1.0, n, nfill, n, \
                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]) - ops.T0

def filled_rolldash(n):
    dn=0.01
    return (logT3(1.0, n+dn, nfill, n+dn, \
                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]) - \
            logT3(1.0, n-dn, nfill, n-dn, \
                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]))/(2.0*dn)

# Single slab of refractive index n observed using converging beam
# note this is T not logT and is used only for refinement
def convslab(n):
    return Tconv(1.0, n, ops.kD[0], ops.freqk) - ops.T0

def convslabdash(n):
    dn=0.01 + 0.001J
    return (Tconv(1.0, n+dn, ops.kD[0], ops.freqk) - \
            Tconv(1.0, n-dn, ops.kD[0], ops.freqk))/(2.0*dn) 

##############################################################################
# Component functions for models included here
##############################################################################
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

def T1(n0, n1, kD):
    "Transmission coefficient of single layer"
    a=4.0*n0*n1/(n0+n1)**2
    b=((n0-n1)/(n0+n1))**2
    return a*exp(-(1.0J)*(n1-n0)*kD)/(1.0 - b*exp(-(2.0J)*n1*kD))

def R1(n0, n1, kD):
    "Reflection coefficient of single layer"
    return ((n0-n1)/(n0+n1))*(1-exp((-2.0J)*n1*kD))/(1-E1(n0, n1, kD))

def T1dash(n0, n1, kD):
    "Derivative of T1"
    return T1(n0, n1, kD)*logT1dash(n0, n1, kD)

def R1dash(n0, n1, kD):
    "Derivative of R1"
    a=-1.0/(n0-n1) - 1.0/(n0+n1)
    b=(2.0J)*kD/(exp((2.0J)*n1*kD) - 1.0)
    return R1(n0, n1, kD)*(a + b + E1dash(n0, n1, kD)/(1.0 - E1(n0, n1,kD)))

def logT1(n0, n1, kD):
    "logT for single layer"
    return log(4.0*n0*n1/(n0+n1)**2) - (1J)*(n1-n0)*kD - log(1.0-E1(n0, n1, kD))

def E1(n0, n1, kD):
    return ((n1-n0)/(n1+n0))**2*exp((-2.0J)*n1*kD)

def logT1dash(n0, n1, kD):
    "Derivative of logT for single layer"
    return 1.0/n1 - 2.0/(n0+n1) - (1J)*kD + E1dash(n0, n1, kD)/(1.0 - E1(n0, n1,kD))

def E1dash(n0, n1, kD):
    return E1(n0, n1, kD)*(2.0/(n1-n0) - 2.0/(n1+n0) + (-2.0J)*kD)

def logT2(n0, n1, n2, kD1, kD2):
    "logT for two layers"
    return logT1(n0, n1, kD1) + logT1(n0, n2, kD2) - \
           log(1 - R1(n0, n1, kD1)*R1(n0, n2, kD2))

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
           log(1 - r1*r0 - r2*r1 - s1*r2*r0)

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

    return logT - log(V)

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
    c=cos(nj*Dj)
    s=sin(nj*Dj)
    d=c + (0.5J)*(nj/n0 + n0/nj)*s
    r=(0.5J)*s*(n0/nj - nj/n0)/d
    t=1.0/d
# note transmission coeff phase is relative to n0 and kD
# s phase which is used to compute reflection coeff is relative to input
    return (r, t*exp((1J)*n0*Dj), t*t-r*r)

def rtsdash(n0, nj, Dj):
    "Derivative of r and t for layer j"
    c=cos(nj*Dj)
    s=sin(nj*Dj)
    (r,t,ss)=rts(n0, nj, Dj)
    d=c + (0.5J)*(nj/n0 + n0/nj)*s
    ddash= (-Dj + (0.5J)*(1.0/n0 - n0/nj**2))*s + (0.5J)*Dj*(nj/n0 + n0/nj)*c
    rdash= r*((nj**2+n0**2)/(nj*(nj**2-n0**2)) + Dj*c/s - ddash/d)
    tdash= -t*ddash/d
#    sdash= 2.0*(tdash*t - rdash*r)
# sdash written this way so it is independent of the phase of t
    sdash= 2.0*((ss+r*r)*tdash/t - rdash*r)
    return (rdash, tdash, sdash)

def one_overTdashdash(n0, n1, kD):
    "(1/T')', used for extending impulse response"
    dn=0.01
    return (1.0/T1dash(n0,n1+dn,kD)-1.0/T1dash(n0,n1-dn,kD))/ \
           (2.0*dn)

############################################################################
# Functions for creating simulated data
############################################################################

def simulate_reference(L, deltat):
    toff=1.0e-11
    twidth=8.0e-13
    tdecay=1.0e-12
    scale=1.0e12
    x=zeros(L, float)
    for t in range(0,L):
        T=t*deltat - toff
        x[t]=-scale*T* exp(-(T/twidth)**2 - T/tdecay)

    return x

def simulate_parallel(x, layers, deltat, noise_level):
    L=len(x)
# convolution done with 4*L points
    M=2*L
    N=4*L
# wave number for frequency point k is k*dk
    deltaf=1.0/(N*deltat)
    dk=2*pi*deltaf/2.9979e8
    m=len(layers)
    indices=    array([l[0] for l in layers])
    thicknesses=array([l[1] for l in layers])
    T=zeros(N, complex)
    for k in range(0,M+1):
        kD=k*dk*thicknesses
        nkD=array([(indices[l],kD[l]) for l in range(0,m)])
        T[k]=RTm(m, 1.0, nkD)[1]
    for k in range(1,M): T[N-k]=conjugate(T[k])
    z=zeros(N, float)
    z[0:L]=x[0:L]
    X=fft(z)/N
    Y=T*X
    y=N*ifft(Y).real + noise_level*random.randn(N)
    return (T,y)

def simulate_converging(x, n0, n1, thickness, deltat, local_theta0, noise_level, shape):
    "wrapper for generating simulated converging data"
    L=len(x)
# use 2*data length to do convolution, note dk is NOT value used in processing
    M=2*L
    N=4*L
    deltaf=1.0/(N*deltat)
# wave number for frequency point k is k*dk
    dk=2*pi*deltaf/2.9979e8
    T=converging(n0, n1, thickness, local_theta0, dk, N, shape)
    z=zeros(N, float)
    z[0:L]=x[0:L]
    X=fft(z)/N
    Y=T*X
    y=N*ifft(Y).real + noise_level*random.randn(N)
    return (T,y)

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
xi=array(xi16+list(-array(xi16)))
wi=array(wi16+wi16)

# use 32 point gauss formula, slower but more accurate
#xi=array(xi32+list(-array(xi32)))
#wi=array(wi32+wi32)

def converging(n0, n1, thickness, theta0, dk, N, shape):
    "used in simulating converging data"
    M=N/2
    local_kscale=dk
    local_shape=shape
# theta_upper_limit is upper limit of theta used in numerical integration
    theta_upper_limit=2.5*theta0
#set up arrays of thetas to be used in integration
    theta=0.5*theta_upper_limit*(1.0 + xi)
    weight=0.5*wi*theta_upper_limit \
        *angular_cutoff(theta, theta0, 0.0, local_kscale, local_shape)
# compute reflection coeffs for these angles, these are independent of freq
    theta_i=theta
    theta_t=arcsin(n0*sin(theta_i)/n1)
    (R1_para, R1_perp, T1_para, T1_perp)=fresnel_oblique(n0, n1, theta_i, theta_t)
    (R2_para, R2_perp, T2_para, T2_perp)=fresnel_oblique(n1, n0, theta_t, theta_i)
#
    T_para = T1_para * T2_para
    R_para = R2_para * R2_para
    T_perp = T1_perp * T2_perp
    R_perp = R2_perp * R2_perp
    s=sin(theta)
    c=cos(theta)
#
    TF=zeros(N,complex)
    i=0
    for k in range(0,M+1):
        kD=k*dk*thickness
# P3, P4 are arrays over theta for this freq
        P4=exp(-1.0J*kD*n1*cos(theta_t))
        P3=exp(-1.0J*kD*(n1*cos(theta_t) - n0*cos(theta_i)))
        F_para=T_para*P3/(1.0 - P4*P4*R_para)
        F_perp=T_perp*P3/(1.0 - P4*P4*R_perp)
        TF[i]=sum((c*c*F_para + F_perp)*s*weight)
        i=i+1
    for i in range(1,M): TF[N-i]=conjugate(TF[i])
    norm=sum((1.0 + c*c)*s*weight)
    return TF/norm

def fresnel_oblique(ni, nt, theta_i, theta_t):
    "Fresnel reflection and transmission coefficients for oblique incidence"
    ci=cos(theta_i)
    ct=cos(theta_t)
    d_perp=ni*ci + nt*ct
    d_para=ni*ct + nt*ci
    r_perp=(ni*ci-nt*ct)/d_perp
    r_para=(nt*ci-ni*ct)/d_para
    t_perp=2.0*ni*ci/d_perp
    t_para=2.0*ni*ci/d_para
    return (r_para, r_perp, t_para, t_perp)

def angular_cutoff(theta, theta0, freqk, kscale, shape):
    "model for cutting off angular spectrum at theta0"
    if shape=="gaussian":
        return exp(-(theta/theta0)**2)
    if shape=="square":
        return 1.0/(1.0 + ((theta/theta0)**6)*(1.0 + (freqk*kscale*theta0)**6))

def set_up_converging_fit(theta0, freqk):
    "creates arrays theta and weight for fitting function Tconv"
    global theta , weight, sintheta, costheta, fit_norm
    local_kscale=0
    local_shape="gaussian"
    theta_upper_limit=2.5*theta0
    theta=0.5*theta_upper_limit*(1.0 + xi)
    weight=0.5*wi*theta_upper_limit* \
      angular_cutoff(theta, theta0, freqk, local_kscale, local_shape)
    sintheta=sin(theta)
    costheta=cos(theta)
    fit_norm=sum((1.0 + costheta*costheta)*sintheta*weight)

def set_angular_spectrum_params(theta0, kscale, shape):
    "used to get parameters into module 'models'"
    global local_theta0, local_kscale, local_shape
    local_theta0=theta0
    local_kscale=kscale
    local_shape=shape

def Tconv(n0, n1, kD, freqk):
    "fitting function for converging beam"
#    set_up_converging_fit(local_theta0, freqk)
    theta_t=arcsin(sintheta/n1)
    costheta_t=cos(theta_t)
    sintheta_t=sin(theta_t)
    sinplus = sintheta*costheta_t + costheta*sintheta_t
    cosplus = costheta*costheta_t - sintheta*sintheta_t
    sinminus= sintheta*costheta_t - costheta*sintheta_t
    cosminus= costheta*costheta_t + sintheta*sintheta_t
    r_perpsq=(sinminus/sinplus)**2
    r_parasq=r_perpsq*(cosplus/cosminus)**2
    P4=exp(-1.0J*kD*n1*costheta_t)
    P3=exp(-1.0J*kD*(n1*costheta_t - n0*costheta))
    t_perp12=4.0*costheta*sintheta_t*costheta_t*sintheta/sinplus**2
    t_para12=t_perp12/cosminus**2
    F_para=t_para12*P3/(1.0 - P4*P4*r_parasq)
    F_perp=t_perp12*P3/(1.0 - P4*P4*r_perpsq)
    return sum((costheta*costheta*F_para + F_perp)*sintheta*weight)/fit_norm


