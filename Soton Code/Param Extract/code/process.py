# version 3.0 25/10/12
import sys
from numpy import *
import ops
import models
from etc import conjgrad, noexit, error_band
from numpy.fft import fft, ifft
from physics import graph, line, clear_graph, pause, postscript

# Threelayer sandwich, fitting middle layer
#def GBsand(nfill):
#    cell=array([(nbread,ops.kD[0]), (nfill,ops.kD[1]), (nbread,ops.kD[2])])
#    return models.RTm(3, 1.0, cell)[1] - ops.T0

#def GBsanddash(nfill):
#    cell=array([(nbread,ops.kD[0]), (nfill,ops.kD[1]), (nbread,ops.kD[2])])
#    return models.RTmdash(3, 1.0, cell, 1)

# Filled roll, fitting outer layers
#def filled_roll(nbread):
#    cell=array([(nbread,ops.kD[0]), (nfill,ops.kD[1]), (nbread,ops.kD[2])])
#    return logT3(1.0, nbread, nfill, nbread, \
#                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]) - ops.T0

#def filled_rolldash(nbread):
#    cell=array([(nbread,ops.kD[0]), (nfill,ops.kD[1]), (nbread,ops.kD[2])])
#    dn=0.01
#    return (logT3(1.0, nbread+dn, nfill, nbread+dn, \
#                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]) - \
#            logT3(1.0, nbread-dn, nfill, nbread-dn, \
#                 thickscale*ops.kD[0], ops.kD[1], thickscale*ops.kD[2]))/(2.0*dn)

# This file contains experimental code for extending the impulse response
#execfile("newsmooth.py")
##########################################################################
# MAIN PROGRAM
############################################################################
clear_graph()
# clear any transfer function used to generate synthetic data
T=[]
#--------------------------------------------------------------------------
# EXECUTE THE CODE TO LOAD DATA
#--------------------------------------------------------------------------
#
# the next line lets one use the command line to specify data
#execfile(sys.argv[1])
#
execfile("example.py")

#
# At this point the data is in arrays x and y of length L and error in y in dy
# M is the length over which impulse response is computed
# N is the zero padded length used for convolutions and display
L=len(x); M=2*L; N=4*L
deltaf=1.0/(N*deltat)
# wave number for frequency point k is k*dk
dk=2*pi*deltaf/2.9979e8
# make sure thickness is an array
thickness=array(thickness)
#-------------------------------------------------------------------------
# Ignore sections of data
#-------------------------------------------------------------------------
if len(ignore)>0: dy=ops.ignore_data(dy, ignore)
if display & 1:
    print "Raw data"
    graph(x[0:L]); graph(y[0:L]); pause(); clear_graph()
#--------------------------------------------------------------------------
# Zero fill data and fft
#--------------------------------------------------------------------------
# The data is ALWAYS extended with an equal length of zeros
x=hstack((x, zeros(L, float)))
y=hstack((y, zeros(L, float)))
# X and Y are fft of data of length M
X=fft(x)/M; Y=fft(y)/M; Xconj=conjugate(X)
#--------------------------------------------------------------------------
# Various conversions to impulse response
#--------------------------------------------------------------------------
# space for selected solution
w=zeros(M,float)

if processing[0]==1:
# direct fft solution with no noise suppression
    display_message="direct fft solution with no noise suppression"
    W=Y/X
# different arrays are used for results so they are availble for comparison
    w1=N*ifft(W).real
    w=copy(w1)
# compute error bars
    dY=fft(dy)/M
    dW=dY/X
    dw=M*ifft(dW).real

if processing[0]==2:
# conjugate gradient method with smoothing
    display_message="conjugate gradient fit to impulse response"
# xx and yy are data expanded to N=2M points
# ee set to ignore data which is beyond recorded end and zero padded to N
    xx=hstack((x, zeros(M, float)))
    yy=hstack((y, zeros(M, float)))
    ee=hstack((1.0/dy**2, zeros(L+M, float)))
# XX and YY are ffts of zero padded data of length N
    XX=fft(xx)/N; YY=fft(yy)/N; XXconj=conjugate(XX)
    ZZ=fft(yy*ee)/N
    exy=(N*N*ifft(ZZ*XXconj).real)[0:M]
    w2=ops.impulse_fit(XX, XXconj, yy, ee, ops.applyEXX, exy, zeros(M, float))
# At this point w2 is causal impulse length M
    w=copy(w2)

# compute error bars
    EE=fft(ee)/N
    ZZ=fft(xx**2)/N
    z=N*ifft(N*EE*conjugate(ZZ)).real
    dw=zeros(M, float)
    dw[0:L]=1.0/sqrt(abs(z[0:L]))
    rmsw=sqrt(sum(w**2)/N)
# fix to prevent large errors where w is not measured
    dw=dw*rmsw**2/(rmsw**2 + dw**2)

#--------------------------------------------------------------------------
# Display impulse response
#--------------------------------------------------------------------------
# At this point we have a fitted impulse response of length M.
if display & 2:
    print display_message
    if display_errors: error_band(w, dw, False)
    else: graph(w)
    if display_theory: graph((ifft(T).real)[0:M])
    pause(); clear_graph()

#--------------------------------------------------------------------------
# Look at Residuals
#--------------------------------------------------------------------------
# look at residuals in time
if display & 4:
# expand w to length N for convolution using fft
    ww=hstack((w, zeros(M, float)))
    WW=fft(ww)/N
    ypred=N*N*ifft(WW*XX).real[0:L]
    diff=yy[0:L] - ypred[0:L]
    print "Data(red), Residuals in time(yellow) and assumed error(green)"
    graph(yy[0:L]); graph(10*abs(diff[0:L])); graph(10*abs(dy[0:L]))
    pause(); clear_graph()

# look at residuals in freq
if display & 8:
    DIFF=YY-N*WW*XX
    print "Data(red), Residuals in frequency(yellow)"
    graph(abs(YY)[0:M]); graph(10*abs(DIFF)[0:M])
    pause(); clear_graph()

#--------------------------------------------------------------------------
# Generate definitive transfer function for conversion to ri
#--------------------------------------------------------------------------
z=hstack((w, zeros(M, float)))
# Note this N*(fft(z)/N) so that transparent has |W|=1
W=fft(z)
# transform error to frequency domain
z=hstack((dw, zeros(M, float)))
dW=fft(z)
# compute sig which is used to compute range of good data
z=hstack((x, zeros(M, float)))
Z=fft(z)/N
sig=abs(Z[0:M])

if display & 16:
    print "Transfer function"
    if display_errors:
        error_band(abs(W[0:M]), abs(dW[0:M]), False)
    else: graph(abs(W)[0:M])
    if display_theory: graph(abs(T))
    pause(); clear_graph()
#---------------------------------------------------------------------------
# Look at reference spectrum to estimate range of good signal.
#---------------------------------------------------------------------------
#sig=abs(X[0:L])
kstart=argmax(sig[20:])
kend1=kstart; kend2=kstart
while sig[kend1]>0.1*sig[kstart]: kend1-=1
while sig[kend2]>0.1*sig[kstart]: kend2+=1
#print kend1, kend2
kend1=max(kend1, 0)
kstart=(kend1+kend2)/2
weight=ops.good_signal(N, kstart, 0.3*(kend2-kend1))
# kend is used for upper freq limit of displays
kend=2*kend2

# botch to set reasonable kend
#kstart=250
#kend=1024

#---------------------------------------------------------------------------
# Take log and unwrap phase
#---------------------------------------------------------------------------
# log-amplitude and phase of W
logW=log(W[0:M])
logWfixed=ops.fixphase(logW, kstart, weight)
if display & 32:
    print "Log amplitude and phase of transfer function and weight"
    graph(logW.real[0:kend]+2); graph(0.005*logWfixed.imag[0:kend])
    graph(0.1*weight[0:kend]); pause(); clear_graph()
#--------------------------------------------------------------------------
# Fit model to get refractive index
#--------------------------------------------------------------------------
# use initial guess at thickness
models.set_thickscale(1.0)
if processing[1]==1:
# simple slab
    model=models.slab1; deriv=models.slab1dash
    (k_lo, k_hi, nraw)=ops.convert_to_ri(kstart, thickness, model, deriv, \
                                     logWfixed, dk)
    display_message = "fitting parallel slab"
    print "convergence failed", k_lo, k_hi

if processing[1]==2:
# GBsandwich, fit ri of filling
    models.set_nbread(nquartz)
    model=models.GBsand; deriv=models.GBsanddash
    (k_lo, k_hi, nraw)=ops.convert_to_ri(kstart, thickness, model, deriv, \
                                     logWfixed, dk)
#    display_message = "fitting GB sandwich"
#    print "convergence failed", k_lo, k_hi

if processing[1]==3:
# Filled Roll, fit ri of bread
    models.set_nfill(1.0)
    model=models.filled_roll; deriv=models.filled_rolldash
    (k_lo, k_hi, nraw)=ops.convert_to_ri(kstart, thickness, model, deriv, \
                                     logWfixed, dk)
    display_message =  "fitting filled roll"
    print "convergence failed", k_lo, k_hi

#----------------------------------------------------------------------------
# Display refractive index
#---------------------------------------------------------------------------
if processing[1]<10 and (display & 64):
# magnitude of spectrum and fitted real and imag parts of refractive index
# imag part with flipped sign to fit on graph and scaled up 10x
    print display_message
    print "Amplitude of transfer function, real ri, imag, ri"
# set W[0] to fix vertical scale
    P=0; Q=kend; W[0]=4.0
    graph(abs(W[P:Q])); graph(nraw[P:Q].real); graph(-10*nraw[P:Q].imag)
    pause(); clear_graph()
    1/0
#----------------------------------------------------------------------------
# Processing for converging beam
#---------------------------------------------------------------------------
if processing[1]==10:
# slab in converging beam
# first fit parallel beam
    model=models.slab1; deriv=models.slab1dash
    (k_lo, k_hi, nraw)=ops.convert_to_ri(kstart, thickness, model, deriv, \
                                     logWfixed, dk)
    print "planar convergence failed", k_lo, k_hi

# set up arrays for fitting converging beam
    models.set_up_converging_fit(theta_cutoff, 1e10)
# fit converging beam
    print "fitting converging beam, please wait"
    model=models.convslab; deriv=models.convslabdash
    (k_lo, k_hi, nrefined)=ops.refine_ri(kstart, thickness, nraw, model, \
                                     deriv, W, dk)
    print "refinement convergence failed", k_lo, k_hi
# display results for converging beam
    if display & 64:
        print "Amplitude of transfer function, real ri, imag, ri"
        print "raw fit using parallel beam"
# set W[0] to fix vertical scale
        P=0; Q=kend; W[0]=4.0
        graph(abs(W[P:Q])); graph(nraw[P:Q].real); graph(-10*nraw[P:Q].imag)
        pause()

        print "refined fit using converging beam"
        graph(nrefined[P:Q].real); graph(-10*nrefined[P:Q].imag)
        pause(); clear_graph()
# replace nraw by nrefined for consistency with later processing
    nraw=nrefined
#----------------------------------------------------------------------------
# Further processing
#---------------------------------------------------------------------------

if processing[2]==1:
# additional calculations to get error bars on ri
    logW1=log(W[0:M] + dW[0:M])
    logW2=log(W[0:M] - dW[0:M])
    logW1fixed=ops.fixphase(logW1, kstart, weight)
    logW2fixed=ops.fixphase(logW2, kstart, weight)
    (k_lo, k_hi, nraw1)=ops.convert_to_ri(kstart, thickness, model, deriv, \
                                     logW1fixed, dk)
    (k_lo, k_hi, nraw2)=ops.convert_to_ri(kstart, thickness, model, deriv, \
                                     logW2fixed, dk)
    dn=nraw2-nraw1
#take absolute values of real and imag parts of dn
    dn=abs(dn.real) + (1J)*abs(dn.imag)
    if display & 128:
# graph to set scale
        print "Amplitude of transfer function, real ri, imag, ri"
        P=0; Q=kend; W[0]=4.0
        graph(abs(W[P:Q]))
        error_band(nraw[P:Q].real, dn[P:Q].real, False)
        error_band(-10*nraw[P:Q].imag, -10*dn[P:Q].imag, True)
        pause()

if processing[2]==2:
# Fit ri to raw data
    (k_lo, k_hi, nrefined)=ops.refine_ri(kstart, thickness, nraw, \
                           models.rawslab1, models.rawslab1dash, W, dk)
    print "fitting ri to raw data"
    print "refinement convergence failed", k_lo, k_hi
    if display & 128:
        print display_message
        print "Amplitude of transfer function(red)"
        print "original ri: real(yellow), imag(green)"
        print "refined ri: real(blue), imag(purple)"
# set W[0] to fix vertical scale
        P=0; Q=kend; W[0]=4.0
        graph(abs(W[P:Q]))
        graph(nraw[P:Q].real); graph(-10*nraw[P:Q].imag)
        graph(nrefined[P:Q].real); graph(-10*nrefined[P:Q].imag)
        pause(); clear_graph()

#if processing[2]==3:
# additional processing to reduce oscillations in ri
#    print "This option is still experimental"
#    pause()
#    Lambda=1e-7
#    set_up_extend()
# check chisq
#    z=hstack((w[0:M], zeros(M, float)))
#    Z=fft(z)/N
#    r=yy - N*N*ifft(Z*XX).real
#    print "original chisq", sum(ee*r*r)
# loop
#    itnum=0
#    while(True):
#        d1=compute_d1()
# solve equations
#        dw=conjgrad(applyd2, d1, zeros(M, float), noexit, 3, True)[0]
# construct modified impulse response
#        neww=neww+dw
# compute new chisq
#        z[0:M]=w[0:M] + dw
#        Z=fft(z)/N
#        r=yy - N*N*ifft(Z*XX).real
#        print "new chisq", sum(ee*r*r)
# compute new transfer function
#        z=hstack((neww[0:M], zeros(M, float)))
#        newW=fft(z)[0:M]
#        newlogW=log(newW[0:M])
#        newlogWfixed=ops.fixphase(newlogW, kstart, weight)
#        model=models.slab1; deriv=models.slab1dash
# fit ri to new transfer function
#        (k_lo, k_hi, n)=ops.convert_to_ri(kstart, thickness, model, deriv, \
#                                          newlogWfixed, dk)
#        print "convergence failed", k_lo, k_hi
# compute smoothness measure before and after
#        smooth1=0; smooth2=0
#        for k in Krange:
#            smooth1+=gamma[k]*abs(nraw[k+1] - nraw[k])**2
#            smooth2+=gamma[k]*abs(n[k+1] - n[k])**2
#        print "roughness", smooth1, smooth2
#        print "dwsq", sum(dw**2)
#        itnum+=1
# emergency exit if loops forever
#        if itnum>10: break
# normal exit from loop
#        if sum(dw**2)<=0.002*sum(w**2): break
#
#    if display & 128:
#        print "Original ri values yellow, green"
# set W[0] to fix vertical scale
#        P=0; Q=kend; W[0]=4.0
#        graph(abs(W[P:Q])); graph(nraw[P:Q].real); graph(-10*nraw[P:Q].imag)
#        print "Improved ri values blue, purple"
#        graph(n[P:Q].real); graph(-10*n[P:Q].imag)
#        pause(); clear_graph()
#        print "Impulse responses"
#        if display_theory:
#            print "theory(red), extended(yellow), truncated(green)"
#            graph((ifft(T).real)[0:M]); graph(neww[0:M]); graph(w[0:L]);
#        else:
#            print "extended(yellow), truncated(green)"
#            graph(neww[0:M]); graph(neww[0:M]); graph(w[0:L]);
#        pause()

if processing[2]==4:
    print "This option is VERY experimental"
# set up various k ranges
    Krange=range(k_lo+2, k_hi-2)
    K_all=slice(k_lo+2, k_hi-2)
    Kminus=slice(k_lo+1, k_hi-3)
    Kplus=slice(k_lo+3, k_hi-1)
    K_lo=slice(k_lo+1, k_hi-2)
    K_hi=slice(k_lo+2, k_hi-1)
# default gamma includes all frequencies, overwrite below for good signal only
    gamma=ones(N, float)
# alias to keep notation same as notes, weight is already double length
#    gamma=weight
    Gamma=zeros(M, float)
    Gamma[K_all]=gamma[K_all] + gamma[Kminus]
    smooth=0.0
    for k in Krange: smooth+=gamma[k]*abs(nraw[k+1] - nraw[k])**2
    print "original roughness", smooth

    n=nraw
    T=zeros(M, complex)
    Tdash=zeros(M, complex)
# parameter determining importance of fitting data
    mu=50.0
# set dnmax to force entry to loop
    dnmax=1.0

    clear_graph()

# mods for converging beam
    mu=1.0
# use this name to distinguish from dn later
    Dn=0.01 + 0.001J
    for it in range(0,5):
        for k in Krange:
# planar code
            T[k]=models.T1(1.0, n[k], k*dk*thickness[0])
            Tdash[k]=models.T1dash(1.0, n[k], k*dk*thickness[0])
#
# converging code
#
#            T[k]=models.Tconv(1.0, n[k], k*dk*thickness[0], ops.freqk)
#            Tdash[k]=(models.Tconv(1.0, n[k]+Dn, k*dk*thickness[0], ops.freqk) - \
#                      models.Tconv(1.0, n[k]-Dn, k*dk*thickness[0], ops.freqk))/(2.0*Dn) 

        Wdefault=fft(hstack((w, zeros(M, float))))
        Lambda=zeros(L, complex)

        fdash=get_fdash(n)
        g=get_g(n)
        impulsepower=sum(w**2)
        print "constraint equation", sum(g**2), impulsepower
        RHS=get_RHS()
        dn=cmplxconjgrad(LHS, LHS, RHS, zeros(M,complex), 1e-3, 3, False)[0]
        print "mean square change in ri", sqrt(sum(dn*conjugate(dn)).real/M)
        n=n + dn
        dnmax=max(abs(dn))
        print "maximum change in ri", dnmax
        smooth=0.0
        for k in Krange: smooth+=gamma[k]*abs(n[k+1] - n[k])**2
        print "current roughness", smooth

        diff=conjugate(fdash) - applyGdag(Lambda)
        print "minimisation derivative", sum(diff*conjugate(diff)).real, "\n"

        rhs=get_rhs()
        Lambda=conjugate(cmplxconjgrad(lhs, lhs, rhs, Lambda, 1e-3, 3, False)[0])
        if display & 128: graph(abs(dn[0:200]))

    if display & 128:
        print "changes in refractive index with iteration"
        pause()
        clear_graph()

    print "raw refractive index (red and green) improved ri (yellow and blue)"
    P=1; Q=300
    graph(nraw.real[P:Q])
    graph(n.real[P:Q])
    graph(-10*nraw.imag[P:Q])
    graph(-10*n.imag[P:Q])

    pause()
    clear_graph()
    Wext=zeros(N, complex)
# prevents offset on impulse
    Wext[0]=W[0]
    for k in range(1,M):
        Wext[k]=T[k]
        Wext[N-k]=conjugate(T[k])

    print "Reconstructed impulse response (red) and  measured (green)"
    graph(ifft(Wext)[0:M].real)
    graph(w[0:L])
    graph(w[0:L])

# Write results to file
#deltat=1.33333e-13
#deltaf=1.0/(N*deltat)
# file = open('n.txt','w')
# for i in range(0, 250):
#     file.writelines(str(nraw[i].real) + "," + str(nraw[i].imag) + "," + \
#                     str(n[i].real) + "," + str(n[i].imag) + "\n")
# file.close()

# 1/0

# file = open('nimag.txt','w')
# for i in range(0, 1024):
#     file.writelines(str(i*deltaf*1e-12)+","+(str(n[i].imag)+"\n"))
# file.close()

# file = open('amp.txt','w')
# for i in range(0, 1024):
#     file.writelines(str(i*deltaf*1e-12)+","+(str(abs(W[i]))+"\n"))
# file.close()
