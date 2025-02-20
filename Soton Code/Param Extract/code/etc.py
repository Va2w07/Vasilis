from physics import line, graph
from numpy.fft import fft, ifft
#from cmath import *      # note this prevents Ufuncs existing
from models import *

global Z, O

def test_models(N,M,dk):
############################################################################
# Tests on models for consistency                                          #
############################################################################
# Triple sandwich
    n0=1.0; n1=1.5; n2=2.0; n3=1.5; kD1=2.0; kD2=0.5; kD3=2.0
# list of layers, elements numbered 0, 1, 2
    layers=[(n1,kD1), (n2, kD2), (n3, kD3)]

# Tests on Reflection and Transmission coefficients for m layers
    (R,T)=RTm(3, n0, layers)
# derivative wrt middle layer
    dT2=RTmdash(3, n0, layers, 1)

# Coefficients with shifted refractive indices
# approx derivative is (Tp-Tm)/(2*dn)
    dn=0.01
    layers=[(n1,kD1), (n2-dn, kD2), (n3, kD3)]
    (Rm,Tm)=RTm(3, n0, layers)
    layers=[(n1,kD1), (n2+dn, kD2), (n3, kD3)]
    (Rp,Tp)=RTm(3, n0, layers)

# Test derivative for 3 layers
    print("English Sandwich, derivative wrt middle layer")
    print( dT2, (Tp-Tm)/(2*dn))

# Test derivative for 2 layers only first two elements of list used
    print("Danish Sandwich, derivative wrt first layer")
    layers=[(n1,kD1), (n2, kD2), (n3, kD3)]
    (R,T)=RTm(2, n0, layers)
    dT1=RTmdash(2, n0, layers, 0)
    layers=[(n1-dn,kD1), (n2, kD2), (n3, kD3)]
    (Rm,Tm)=RTm(2, n0, layers)
    layers=[(n1+dn,kD1), (n2, kD2), (n3, kD3)]
    (Rp,Tp)=RTm(2, n0, layers)
    print( dT1, (Tp-Tm)/(2*dn))

    print("Danish Sandwich, derivative wrt second layer")
    layers=[(n1,kD1), (n2, kD2), (n3, kD3)]
    (R,T)=RTm(2, n0, layers)
    dT2=RTmdash(2, n0, layers, 1)
    layers=[(n1,kD1), (n2-dn, kD2), (n3, kD3)]
    (Rm,Tm)=RTm(2, n0, layers)
    layers=[(n1,kD1), (n2+dn, kD2), (n3, kD3)]
    (Rp,Tp)=RTm(2, n0, layers)
    print(dT2, (Tp-Tm)/(2*dn))

# Check r and t for each layer against separate code for R1 and logT1
    print("check r and t for each layer")
    layers=[(n1,kD1), (n2, kD2), (n3, kD3)]
    (r0,t0,s)=rts(n0, layers[0][0], layers[0][1])
    (r1,t1,s)=rts(n0, layers[1][0], layers[1][1])
    (r2,t2,s)=rts(n0, layers[2][0], layers[2][1])
    print(r0, R1(n0, n1, kD1))
    print (t0, exp(logT1(n0, n1, kD1)))
    print (r1, R1(n0, n2, kD2))
    print (t1, exp(logT1(n0, n2, kD2)))
    print (r2, R1(n0, n3, kD3))
    print (t2, exp(logT1(n0, n3, kD3)))

# Check transmission of 2 layers against analytic result
    print( "Check combined layer against individual")
    print( T, t0*t1/(1-r0*r1))

# Check logTm against RTm logT1, logT2, logT3 and derivatives
    print(  "Check logT1 for 1 layer")
    (R,T)=RTm(1, n0, layers)
    print( exp(logTm(1, n0, layers)), T)
    print( logTm(1, n0, layers), logT1(n0,n1,kD1))
    print( "Check derivative")
    dn=0.005
    lTp=logT1(n0,n1+dn,kD1)
    lTm=logT1(n0,n1-dn,kD1)
    print (logT1dash(n0,n1,kD1), (lTp-lTm)/(2*dn) )
    print ("Check log2T for 2 layers")
    (R,T)=RTm(2, n0, layers)
    print (exp(logTm(2, n0, layers)), T)
    print( logTm(2, n0, layers), logT2(n0,n1,n2,kD1,kD2))
    print( "Check derivative")
    lTp=logT2(n0,n1,n2+dn,kD1,kD2)
    lTm=logT2(n0,n1,n2-dn,kD1,kD2)
    print( logT2dash(n0,n1,n2,kD1,kD2), (lTp-lTm)/(2*dn) )
    print ("Check logT3 for 3 layers")
    (R,T)=RTm(3, n0, layers)
    print (exp(logTm(3, n0, layers)), T)
    print (logTm(3, n0, layers), logT3(n0,n1,n2,n3,kD1,kD2,kD3))
    print( "Check derivative")
    lTp=logT3(n0,n1,n2+dn,n3,kD1,kD2,kD3)
    lTm=logT3(n0,n1,n2-dn,n3,kD1,kD2,kD3)
    print( logT3dash(n0,n1,n2,n3,kD1,kD2,kD3), (lTp-lTm)/(2*dn) )

# Check layers of same refractive index against thicker layers
    n1=2.0
    kD=1.5
    layers=[(n1,kD),(n1,kD),(n1,kD)]
    RTA=RTm(3, 1.0, layers)
    layers=[(n1,3*kD)]
    RTB=RTm(1, 1.0, layers)
    print ("Reflection coefficient for 3 layers and 1 layer 3 time thicker")
    print (RTA[0], RTB[0])
    print ("Transmission coefficient for 3 layers and 1 layer 3 time thicker")
    print (RTA[1], RTB[1])
# Check logT3 against logT1
    TA=logT3(1.0,2.0,2.0,2.0,kD,kD,kD) 
    TB=logT1(1.0,2.0,3*kD)
    print( "Check logT3 against logT1 for same index")
    print( TA, TB)

# Consistency check on generating and analysis functions for GBsand
    layers=[(2.0, 1.0e-3), (3.0-0.02J, 0.5e-3), (2.0, 1.0e-3)]
    nbread=layers[0][0]
    nfill= layers[1][0]
    thickness=array([layers[0][1], layers[1][1], layers[2][1]])
    indices=    array([l[0] for l in layers])
    thicknesses=array([l[1] for l in layers])
    m=len(layers)
    T1=zeros(N, complex)
    for k in range(0,M):
        kD=k*dk*thicknesses
        nkD=array([(indices[l],kD[l]) for l in range(0,m)])
        T1[k]=RTm(m, 1.0, nkD)[1]

    T2=zeros(N, complex)
    for k in range(0,N):
        kD=k*dk*thickness
        T2[k]=exp(logT3(1.0,nbread,nfill,nbread,kD[0],kD[1],kD[2]))

    graph(T1[0:M].real)
    graph(T2[0:M].real)
    graph(T1[0:M].imag)
    graph(T2[0:M].imag)

############################################################################
# Test on converging transmission
############################################################################

def test_converging(N, M, dk):
    global theta_cutoff, n1, n2, kD
    theta_cutoff=0.5
    theta_max=2.0*theta_cutoff
    norm0=gaussint(0.0, theta_max, 1e-6, f0)
    n1=1.0
    n2=2.0
    thickness=1.0e-3

# check value of transmission at selected frequencies against integral
    print ("Check transmission coefficient at selected frequencies")
    print ("Frequency point 100")
    freq_num=100
    kD=freq_num*dk*thickness
    print (gaussint(0.0, theta_max, 1e-6, f)/norm0)
    print (converging(n1, n2, 1.0e-3, theta_cutoff, dk, N)[freq_num])

    print( "Frequency point 200")
    freq_num=200
    kD=freq_num*dk*thickness
    print( gaussint(0.0, theta_max, 1e-6, f)/norm0)
    print (converging(n1, n2, 1.0e-3, theta_cutoff, dk, N)[freq_num])

# Check transmission for very small theta_cutoff
    TC=zeros(N,complex)
    TP=zeros(N,complex)
    theta_cutoff=0.02
    for k in range(0,M):
        kD=k*dk*thickness
        TP[k]=T1(n1, n2, kD)
    for k in range(1,M): TP[N-k]=conjugate(TP[k])

    TC=converging(n1, n2, 1.0e-3, theta_cutoff, dk, N)
    print ("Graph checks transmission coefficient for very small cutoff")
    graph(abs(TP[0:M]))
    graph(abs(TC[0:M]))

# Test consistency of functions for generating and analysing converging data
    n1=2.0
    thickness=2.0e-3
    theta_cutoff=0.7
    TA=converging(1.0, n1, thickness, theta_cutoff, dk, N)
    set_up_converging_fit(theta_cutoff)
    print ("Check creation and fitting of converging data")
    print ("frequency point 100")
    freq_num=100
    kD=freq_num*dk*thickness
    TB=Tconv(1.0, n1, kD)
    print( TA[freq_num], TB)
    print( "frequency point 200")
    freq_num=200
    kD=freq_num*dk*thickness
    TB=Tconv(1.0, n1, kD)
    print( TA[freq_num], TB)

def f0(theta):
    s=sin(theta)
    c=cos(theta)
    return angular_cutoff(theta, theta_cutoff)*(1+c*c)*s

def f(theta):
# first interface
    theta1=theta
    theta2=arcsin(n1*sin(theta1)/n2)
    (R1_para, R1_perp, T1_para, T1_perp)=fresnel_oblique(n1, n2, theta1, theta2)
# second interface
    (R2_para, R2_perp, T2_para, T2_perp)=fresnel_oblique(n2, n1, theta2, theta1)
# propagation
    P2=exp(-1.0J*kD*n2/cos(theta2))
    P1=exp(-1.0J*kD*n1/cos(theta1))
    F_para=T2_para*(P2/P1)*T1_para/(1.0 - P2*P2*R2_para*R2_para)
    F_perp=T2_perp*(P2/P1)*T1_perp/(1.0 - P2*P2*R2_perp*R2_perp)
    s=sin(theta)
    c=cos(theta)
    return angular_cutoff(theta, theta_cutoff)*(c*c*F_para + F_perp)*s

############################################################################
# Various functions for computing transmission through systems
#############################################################################
def set_up_frequencies(N, fmax):
    "Returns array of N frequencies from -fmax to +fmax, in units of fref"
    "The frequencies are in FFT order (negative after positive)"
    "Also sets global variables O and Z for use elsewhere"
    global O, Z
    Z=zeros(N,float)
    O=ones(N,float)
    f=zeros(N,float)
    for q in range(0,N/2): f[q]=1+(2*q*fmax)/N
    for q in range(N/2,N): f[q]=1+(2*(N-q)*fmax)/N
    return f

def wavenumber(f, nsq):
    "Computes wavenumbers for frequencies in list f using values of the"
    "squares of the corresponding refractive indices in nsq"
    "The sign of the square root is chosen to make the imaginary part <0"
    "Note k in program is ck/Fref in theory"
# Note can't take complex functions of arrays, hence loop)
    klist=[]
    for w,n in zip(f, nsq):
        k=2*pi*w*sqrt(n)
        if k.imag > 0: k=-k
        klist.append(k)

    return array(klist)

#def gaussian_pulse(N, width, centre):
#    "Generates array containing gaussian pulse of total length"
#    " N points with given width and centre" 
#    return array([exp(-((n-centre)/width)**2) for n in range(0,N)], Complex)

def lorentz_susceptibility(omega, A, wr, wp):
    "Computes lorentz susceptibility of medium as function of freq"
    return array([wp*wp/(-w*w + (2J)*A*w + wr*wr) for w in omega])

def propagate(klist, D, direction):
    "Computes propagator component"
    forward=expikD(klist, D, -direction)
    return (forward, Z, Z, forward)

def expikD(klist, D, sign):
    "Computes exp(+-ikD)"
    return array([exp(sign*(1J)*k*D)  for k in klist])

def fresnel(k1, k2):
    "Computes transmission and reflection component for interface"
    "Wave goes from medium k1 to medium k2"
    R=(k1-k2)/(k1+k2)
    return (1+R, -R, R, 1-R)

def slab(k_inside, k_outside, D):
    "Computes transmission and reflection coefficients for slab of"
    "thickness D of medium with wavenumber k_inside, embedded in"
    "medium with wavenumber k_outside"
    T=[]; R=[]
    for k_in, k_out in zip(k_inside, k_outside):
        s=sin(k_in*D)
        c=cos(k_in*D)
        d=4*k_in*k_out*c + (2J)*(k_in*k_in + k_out*k_out)*s
        T.append(4*k_in*k_out/d)
        R.append((2J)*s*(k_out**2 - k_in**2)/d)

    T=array(T); R=array(R)
    return (T, R, R, T) 

def mirror(omega):
    "Computes perfect mirror component"
    return (Z, Z, -O, Z)

#
# A 4-port object is described by the two outputs as a function of the two
# inputs.  It is convenient to split this into a linear and non-linear part
#
# O1 = L11 * I1 + L12 * I2
# O2 = L21 * I1 + L22 * I2
#
# where L11, L12, L21, L22 are arrays describing the frequency dependence.
# As a python object these are packed as a tuple: (L11, L12, L21, L22)
#
def concat(first, second):
    "Adds 4-port non-linear device to another a-port non-linear device"
    " returns 4-port non-linear device"
    (L11, L12, L21, L22)=first
    (l11, l12, l21, l22)=second
    N=len(L11)
    F=1/(1 - L12*l21)
    T11=l11*F*L11
    T12=l11*F*L12*l22 + l12
    T21=L21 + L22*F*l21*L11
    T22=L22*F*l22
    return (T11, T12, T21, T22)


def transit(ew1, ew2, component):
    "Computes E field (in freq) after transiting component"
    "Inputs ew1 and ew2 are fields entering from left and right"
    "Returns fields at left and right on exit"
    (L11, L12, L21, L22, N12)=component
    return(L11*ew1+L12*ew2, L21*ew1+L22*ew2)


def gaussint(a, b, delta, f):
    "Integrates f(x) between a and b with accuracy delta"
    n=8
    f1=I(n, a, b, f)
    f2=-f1
    while abs(f2-f1) > delta*abs(f1):
        if (n==16384):
            print ("Can't achieve accuracy requested")
            return 0

        f2=f1
        n=2*n
        f1=I(n, a, b, f)

    return f1

#---------------------------------------------------------------------------

def I(n, a, b, f):
    "Divides range a to b into n panels"
    h=(b-a)/n
    sum=0.0
    for i in range(0,n):
        sum+=tenpoint(a+i*h, a+(i+1)*h, f)
    return sum

def tenpoint(a, b, f):
    "Evaluates integral from a to b by 10 point Gauss formula"
    u=0.5*(a+b); v=0.5*(b-a)
    s=0.2955242247*(f(u+0.1488743389*v)+f(u-0.1488743389*v))+ \
      0.2692667193*(f(u+0.4333953941*v)+f(u-0.4333953941*v))+ \
      0.2190863625*(f(u+0.6794095682*v)+f(u-0.6794095682*v))+ \
      0.1494513491*(f(u+0.8650633666*v)+f(u-0.8650633666*v))+ \
      0.0666713443*(f(u+0.9739065285*v)+f(u-0.9739065285*v))
    return 0.5*s*(b-a)

############################################################################
# Colour plots of 2D arrays
#############################################################################
def set_up_colours():
    global colour
    colour=[0 for i in range(0,256)]
    for i in range(0, 32):    colour[i] = (mix(i,32), 0, 0)
    for i in range(32, 80):   colour[i] = (255, mix(i-32,48), 0)
    for i in range(80, 128):  colour[i] = (mix(128-i,48), 255, 0)
    for i in range(128, 176): colour[i] = (0, 255, mix(i-128,48))
    for i in range(176, 225): colour[i] = (0, mix(224-i,48), 255)
    for i in range(225, 256): colour[i] = (mix(i-224,32), 0, 255)
    for i in range(0,256): colour[i]="#%02x%02x%02x" % colour[i]

def mix(i, range):
    f=1.0-i/float(range)
    return int(255*(1-f*f*(3-2*f)))

def colourplot(a):
    M=256
    (w,h)=shape(a)
    N=max(w,h)
    amax=a.max(); amin=a.min()
    if N<=M:
        r=N; s=int(M/N); t=1
    if N>M:
        r=M; s=1; t=int(N/M)
    c=Canvas(width=M, height=M, background="grey60")
    c.pack()
    for i in range(0,r):
        for j in range(0,r):
            p=int(i*t); q=int(j*t)
            col=colour[int(255*(a[p,q] - amin)/(amax - amin))]
            c.create_rectangle(i*s, j*s, i*s+s, j*s+s, width=0, fill=col)

set_up_colours()

############################################################################
# Routines from Numerical Recipes, recoded in Python
#############################################################################

def SIGN(a,b):
    if b>0.0: return  abs(a)
    else:     return -abs(a)

def frprmn(p, ftol, func, dfunc):
    ITMAX=400
    EPS=1.0e-10
    fp=func(p); xi=dfunc(p)
    g=-xi
    h=array(g)
    xi=array(g)
    for its in range(1,ITMAX):
        (xmin, fret)=dlinmin(p, xi, func, dfunc)
#        print (fret,xmin)
#        1/0
        if 2.0*abs(fret-fp) <= ftol*(abs(fret) + abs(fp) + EPS): return p
        fp=func(p); xi=dfunc(p)
        gg=sum(g*g)
#        dgg=sum(xi*xi)
        dgg=sum((xi+g)*xi)
        if gg==0.0: return p
        gam=dgg/gg
        g=-copy(xi)
        h=g + gam*h
        xi=copy(h)
    print ("Too many iterations in frprmn")


def dlinmin(p, xi, func, dfunc):
    global pcom, xicom, nrfunc, nrdfun
    TOL=2.0e-4
    nrfunc=func
    nrdfun=dfunc
    pcom=array(p)
    xicom=array(xi)
    ax=0.0
    xx=1.0
#    print xi
    (ax, xx, bx, fa, fx, fb)=mnbrak(ax, xx, f1dim)
#    print (ax,xx,bx,fa,fx,fb)
#    1/0
    (xmin, fret)=dbrent(ax, xx, bx, f1dim, df1dim, TOL)
    xi=xi*xmin
    p+=xi
    return (xmin, fret)

def f1dim(x):
    xt=pcom + x*xicom
    return nrfunc(xt)

def df1dim(x):
    xt=pcom + x*xicom
    df=nrdfun(xt)
    df1=sum(df*xicom)
    return df1

def mnbrak(ax, bx, func):
    GOLD=1.618034
    GLIMIT=100.0
    TINY=1.0e-20
    fa=func(ax); fb=func(bx)
    if fb>fa: (ax, fa, bx, fb)=(bx, fb, ax, fa)
    cx=bx + GOLD*(bx - ax)
    fc=func(cx)
    while fb>fc:
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx - ((bx-cx)*q - (bx-ax)*r)/(2.0*SIGN(max((abs(q-r),TINY)), q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if (bx-u)*(u-cx) > 0.0:
            fu=func(u)
            if fu<fc:
                ax=bx; bx=u; fa=fb; fb=fu
                return (ax, bx, cx, fa, fb, fc)
            else:
                if fu>fb:
                    cx=u; fc=fu
                    return (ax, bx, cx, fa, fb, fc)
            u=cx + GOLD*(cx-bx)
            fu=func(u)
        elif (cx-u)*(u-ulim) > 0.0:
            fu=func(u)
            if fu<fc:
                bx=cx; cx=u; u=cx+GOLD*(cx-bx)
#                (bx, cx, u)=(cx, u, cx+GOLD*(cx-bx))
                (fb, fc, fu)=(fc, fu, func(u))
        elif (u-ulim)*(ulim-cx) >= 0.0:
            u=ulim; fu=func(u)
        else:
            u=cx + GOLD*(cx-bx)
            fu=func(u)
        (ax, bx, cx)=(bx, cx, u)
        (fa, fb, fc)=(fb, fc, fu)

    return (ax,bx,cx,fa,fb,fc)

def dbrent(ax, bx, cx, f, df, tol):
    ITMAX=100
    ZEPS=1e-10
    e=0.0
    a=min((ax, cx))
    b=max((ax, cx))
    v=bx; w=v; x=w
    fx=f(x); fv=fx; fw=fv
    dx=df(x); dv=dx; dw=dv
    for iter in range(1, ITMAX):
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.0*tol1
        if abs(x-xm) <= tol2-0.5*(b-a):
            return (x, fx)
        if abs(e) > tol1:
            d1=2.0*(b-a)
            d2=d1
            if dw!=dx: d1=(w-x)*dx/(dx-dw)
            if dv!=dx: d2=(v-x)*dx/(dx-dv)
            u1=x+d1; u2=x+d2
            ok1=(a-u1)*(u1-b)>0.0 and dx*d1<=0.0
            ok2=(a-u2)*(u2-b)>0.0 and dx*d2<=0.0
            olde=e
            e=d
            if ok1 or ok2:
                if ok1 and ok2:
                    if abs(d1)<abs(d2): d=d1
                    else:               d=d2
                elif ok1: d=d1
                else:     d=d2
                if abs(d) <= abs(0.5*olde):
                    u=x+d
                    if (u-a < tol2) or (b-u < tol2): d=SIGN(tol1, xm-x)
                else:
                    if dx>=0.0: e=a-x
                    else:       e=b-x
                    d=0.5*e
            else:
                if dx>=0.0: e=a-x
                else:       e=b-x
                d=0.5*e
        else:
            if dx>=0.0: e=a-x
            else:       e=b-x
            d=0.5*e
        if abs(d) >= tol1:
            u=x+d
            fu=f(u)
        else:
            u=x+SIGN(tol1, d)
            fu=f(u)
            if fu > fx: return (x, fx)

        du=df(u)
        if fu<=fx:
            if u>=x: a=x
            else: b=x
            (v,fv,dv)=(w,fw,dw)
            (w,fw,dw)=(x,fx,dx)
            (x,fx,dx)=(u,fu,du)
        else:
            if u<x: a=u
            else: b=u
            if fu<=fw or w==x:
                (v,fv,dv)=(w,fw,dw)
                (w,fw,dw)=(u,fu,du)
            else:
                if fu<fv or v==x or v==w: (v,fv,dv)=(u,fu,du)

    print ("Too many iterations in dbrent")
#    return (x,fx)

# Rosenbrock Test function
#def func(x):
#    return 100*(x[1]-x[0]**2)**2 + (1-x[0])**2
#def dfunc(x):
#    return array([-400*x[0]*(x[1]-x[0]**2)+2*x[0]-2,200*(x[1]-x[0]**2)])
#(ax,bx,cx,fa,fb,fc)=mnbrak(0.1, 0.2, func)
#print (ax,bx,cx,fa,fb,fc)
#print dbrent(ax,bx,cx,func,dfunc,1e-10)
#print frprmn(array([0.5,2]), 1e-5, func, dfunc)

def tridiag1(A, b):
    "solves Ax=b where A has array A down diagonal and ones either side"
    n=size(A)
    u=zeros(n, complex)
    gamma=zeros(n, complex)
    u[0]=b[0]/A[0]
    beta=A[0]
    for j in range(1,n):
        gamma[j]=1.0/beta
        beta=A[j] - gamma[j]
        u[j]=(b[j]-u[j-1])/beta

    for j in range(n-2, -1, -1):
        u[j]=u[j] - gamma[j+1]*u[j+1]

    return u

def tridiag(a, b, c, r):
    "solves Au=r where A has array b down diagonal b[0] to b[n-1], a down "
    "the -1 diagonal a[2] to a[n-1] and c the +1 diagonal c[0] to c[n-2]"
    n=size(b)
    u=zeros(n, complex)
    gamma=zeros(n, complex)
    u[0]=r[0]/b[0]
    beta=b[0]
    for j in range(1,n):
        gamma[j]=c[j-1]/beta
        beta=b[j] - a[j-1]*gamma[j]
        u[j]=(r[j] - a[j-1]*u[j-1])/beta

    for j in range(n-2, -1, -1):
        u[j]=u[j] - gamma[j+1]*u[j+1]

    return u

def conjgrad(A, b, x, chilim, itmax, debug):
    "solves Ax=b by conjugate gradients. A is performed by function."
    EPS=1.0e-14
    n=len(b)
    p=zeros(n, float)
    r=zeros(n, float)
    z=zeros(n, float)
    it=0
    r=b-A(x)
    bnorm=sqrt(dot(b,b))
    z=precon(r)
    while it<=itmax:
        it+=1
        bknum=dot(z, r)
        if it==1:
            p=copy(z)
        else:
            bk=bknum/bkden
            p=bk*p + z
        bkden=bknum
        z=A(p)
        akden=dot(z, p)
        ak=bknum/akden
        x=x + ak*p
        r=r - ak*z
        z=precon(r)
        err=sqrt(dot(r,r))/bnorm
        if debug: print (it, err)
        if chilim(x): break
    return (x, it, err)

# noexit is used to force an exit on number of iterations
def noexit(x):
    return False

# precon is included so that conjgrad can be speeded up for large datasets
# using noise supression approximate solution.

#def precon(r):
#    nu=1.0e-7
#    b=zeros(M, float)
#    b[0:L]=r[0:L]
#    B=fft(b)/M
#    B=B/(X*Xconj + nu)
#    return (ifft(B).real)[0:L]

def precon(x):
    return copy(x)

############################################################################
# Extension to graph function to plot error bars
#############################################################################
def error_band(y, dy, add):
# draws y superimposed on error bars from y-dy to y+dy.
# If add==True the previous plot and its axes are kept.
    n=len(y)-1
    if not add:
        clear_graph()
# plot a graph to set axes, then clear it
        graph(y)
        clear_graph()
# draw the error band
    for j in range(0,600):
        x=n*float(j)/600.0
        i=int(x)
        f=x-i
# interpolate values
        yi=y[i]*(1.0-f) + y[i+1]*f
        dyi=dy[i]*(1.0-f) + dy[i+1]*f
# if zero is a reasonable value for y show error bars about this value
        if abs(yi)<abs(dyi): yi=0.0
        y1=yi - dyi
        y2=yi + dyi
        line(x, y1, x, y2, 6)
# replot the graph to superimpose
    graph(y)

