from numpy.fft import fft, ifft, fft2, ifft2
from numpy.linalg import solve

import models
import numpy
import ops
import pylab


def extraction(N, M, dk, deltat, deltaf, reference_file, sample_file, thickness, convthickness, setthetacutoff, processing, kstart, setkmin, paddingfactor, phasetolerance, setangularcentre, display, numRefl, polarisation_configuration, angular_spread_sim, theta):
    #extraction(N, M, dk, deltat, deltaf, reference_file, sample_file, thickness, convthickness, setthetacutoff, processing, kstart, setkmin, paddingfactor, phasetolerance, setangularcentre, display, numRefl, polarisation_configuration):
   
     ##########################################################################
    # MAIN PROGRAM
    ############################################################################
    
    # Set nu factor for noise suppression
    nufac=0.1
    # Set refractive index of "bread" for GB sandwich
    nbread=0.0
    # Set name of file to save simulated data
    simulated_data_file="simulated_data.txt"

    
    # Set file names of reference and sample
    # Set thicknesses of layers in sample in metres.  Note thickness may be a list.
    
    models.set_theta_cutoff(setthetacutoff, setkmin, setangularcentre)
    
    length = int(thickness[0]*1.0e6)
    
    # make sure thickness is an array
    thickness=numpy.array(thickness)
    convthickness=numpy.array(convthickness)
    
    tscale=numpy.zeros(M, complex)
    fscale=numpy.zeros(N, complex)
    for i in range(0, N):
        fscale[i]=i*deltaf*1e-12
    for i in range(0, M):
        tscale[i]=i*deltat*1e12
   
    
    
    start=0
    "mx shifts the sample scan by a number of units (crete data has a shift)"
    mx=ops.findtimeshift(reference_file, sample_file)
    #mx=0
    xx=ops.readdatafile(reference_file, 1.0, paddingfactor)[start:]
    x=numpy.zeros(N,float)
    x[0:M]=xx[0:M]
    if processing[0]==0:
        yy=ops.readdatafile(sample_file, 1.0, paddingfactor)[start:]
        y=numpy.zeros(N,float)    
        y[mx:M+mx]=yy[0:M]
        'use when hard truncating in ops'
        #y[mx:M+mx]=yy[0:M-mx]
          
    #this is for simulate data where there is no need to step
        ##y[0:M]=yy[0:M]
    # Create various types of synthetic data
    
    'this is to test the casuality by moving the pulse backwards or forwards'
    nx=0
    if processing[0]==1:
        layers=[(3.500-0.000J, 0.500e-3)]
        y=models.simulate_parallel2(reference_file, layers, dk, 0.0, paddingfactor, N, M, 'none')
        'put the pulse nx position after actual'
        #y[nx:M+nx]=y[0:M]
        'put the pulse nx position before actual'
        #y[0:M-nx]=y[nx:M]
        #ops.savedatafile(y, simulated_data_file, 1e-9)
        #thickness=numpy.array([layers[0][1]])
    if processing[0]==2:
        layers=[(2.0, 1.0e-3), (3.0-0.02J, 0.1e-3), (2.0, 1.0e-3)]
        y=models.simulate_parallel(reference_file, layers, dk, 1e-5)
        ops.savedatafile(y, simulated_data_file, 1e-9)
        nbread=layers[0][0]
        thickness=numpy.array([layers[0][1], layers[1][1], layers[2][1]])
    if processing[0]==3:
        layers=[(1.000-0.005J, 0.500e-3)]
        theta_cutoff=0.24
        y=models.simulate_converging(reference_file, 1.0, 1.00001-0.000J, 0.500e-3, theta_cutoff, dk, 0.0)
        #ops.savedatafile(y, simulated_data_file, 1.0)
        #thickness=numpy.array([layers[0][1]])
        
        #Converging beam with conical beam profile
    if processing[0]==4:
        #layers=[(1.500-0.005J, 3.000-3)]
        theta_c=0.10+0.01J
        y=models.simulate_converging_conical(reference_file, 1.0, 1.500-0.005J, 3.000e-3, theta_c, dk, 0.0, paddingfactor, N, M)
        #ops.savedatafile(y, simulated_data_file, 1.0)
        #thickness=numpy.array([layers[0][1]])
    
    #---------------------------------------------------------------------------
    # Set noise parameters
    #---------------------------------------------------------------------------
    # noise is empirically raised by factor of 5
    dy=5.0*ops.noise_model_time(y)
    dybar=sum(dy)/N
    dysq=dy*dy
    dysqbar=dybar*dybar
    xsqbar=sum(x*x)/N
    nu=nufac*numpy.sqrt(xsqbar/dysqbar)/N
    e=numpy.zeros(N)
    if processing[2]==4: e[0:M]=1.0/dysq[0:M]
    else:                dysq[0:M]=dysqbar
    
    if display & 2048:
        print "Raw data"
        pylab.figure()
        pylab.title('Raw Data')
        pylab.plot(xx[:])
        pylab.plot(yy[:])
        
    #-------------------------------------------------------------------------
    # Truncate the data for long time delays
    #-------------------------------------------------------------------------
    truncateval=50
    
    'crete_data4'
    
    truncateshift=0
    falloff=30
    
    
    if processing[1]==1:

        Q=N/2
        pylab.figure()
        pylab.subplot(211)
        pylab.plot( y[0:Q], label='y')
        pylab.plot( x[0:Q],label='x')        
        pylab.title('scans')
        pylab.legend()
        pylab.grid(True)
        xold = x[:]
        yold = y[:]
        
        "default"
#        x=ops.truncate(x, numpy.argmax(x), 30)
#        y=ops.truncate(y, numpy.argmax(y), 30)
        "crete_data4"
        #print numpy.argmax(x)
#        x=ops.truncate(x, truncateshift, falloff)
#        y=ops.truncate(y, truncateshift, falloff)
#        x=ops.truncate(x, numpy.argmax(x), falloff)
#        y=ops.truncate(y, numpy.argmax(y), falloff)
        
    
        pylab.subplot(212)
        pylab.plot( y[0:Q]/yold[0:Q], label='y')
        pylab.plot( x[0:Q]/xold[0:Q],label='x')        

        pylab.plot( y[0:Q], label='y')
        pylab.plot( x[0:Q],label='x')        
        pylab.title('truncated scans')
        pylab.legend()
        pylab.grid(True)
        
        
    if display & 1024:
        pylab.figure()
        pylab.title('Time Windowed Data')
        print "Time windowed data"
        pylab.plot(x[0:M])
        pylab.plot(y[0:M])
    
    #--------------------------------------------------------------------------
    # Various conversions to transfer function
    #--------------------------------------------------------------------------
    X=fft(x)/N; Y=fft(y)/N; Xconj=numpy.conjugate(X)
    Xlog=numpy.log(X); Ylog=numpy.log(Y)
    # space for selected solution
    kstart = numpy.argmax(Xlog.real[0:len(Xlog)/2])
    w=numpy.zeros(N,float)
    
    if processing[2]==0:
    # direct fft solution with no noise suppression
        display_message="direct fft solution with no noise suppression"
        W=Y/X
        w0=N*ifft(W).real
        w=numpy.array(w0)
    
    if processing[2]==1:
    # direct fft solution with noise suppression
        display_message="direct fft solution with noise suppression"
        W=Y*Xconj/(X*Xconj + nu*dysqbar)
        w1=N*ifft(W).real
        w=numpy.array(w1)
    
    if processing[2]==2:
    # alternative fft with noise reduction with transparent default
        display_message="alternative fft with noise reduction with transparent default"
    # find peak in impulse response
        W=Y*Xconj/(X*Xconj + nu*dysqbar)
        w2=N*ifft(W).real
        peak=-1e10; place=0
        for t in range(0,N):
            if abs(w2[t])>peak: peak=abs(w2[t]); place=t
    # optimise the position of the peak
        best=1e10; pos=0
        for T in range(place-20,place+20):
            ex=numpy.zeros(N,float); ex[T]=1; EX=fft(ex)
            W=(Y*Xconj+nu*dysqbar*EX)/(X*Xconj + nu*dysqbar)
            Z=Y-W*X
            s=sum(Z*numpy.conjugate(Z))
            Z=W-EX
            s=s+nu*dysqbar*sum(Z*numpy.conjugate(Z))
            if s<best: beat=s; pos=T
    
        ex=numpy.zeros(N,float); ex[pos]=1; EX=fft(ex)
        W=(Y*Xconj+nu*dysqbar*EX)/(X*Xconj + nu*dysqbar)
        w4=N*ifft(W).real
        w=numpy.array(w4)
    
    if processing[2]==3:
    # levinson solution
        display_message="causal solution, constant noise"
        xx=ifft(X*Xconj + nu*dysqbar).real
        xy=N*ifft(Y*Xconj).real
        r=numpy.zeros(2*M-1,float)
        for i in range(0,M):
            r[M-1+i]=xx[i]
            r[M-1-i]=xx[i]
        w3=ops.PYlevinson(r,xy[0:M])
        w[0:M]=w3[0:M]
    
    if processing[2]==4:
    # exact solution with causality and variable noise
        display_message="exact solution with causality and variable noise"
        print "Processing, please wait"
        a=numpy.zeros((N,N), float)
        for t in range(0,M): a[t,t]=e[t]
        for t in range(1,M): a[N-t,N-t]=a[t,t]
        A=ifft2(a); A=A*numpy.outer(X,X); a=fft2(A).real; a=a[0:M,0:M]
        for t in range(0,M): a[t,t]+=nu
        b=ifft(Xconj*fft(e*y)); b=b[0:M].real
        w4=solve(a,b)
        w[0:M]=w4[0:M]
    
    if display & 512:
        pylab.figure()
        pylab.title(display_message)
        print display_message
        pylab.plot(w)
    
    # generate definitive spectrum
    W=fft(w)/N
    #--------------------------------------------------------------------------
    # Look at Residuals
    #--------------------------------------------------------------------------
    # look at residuals in time
    if display & 256:
        ypred=N*ifft(W*X).real
        diff=y[0:M] - ypred[0:M]
        pylab.figure()
        pylab.title("Data, Residuals in time and assumed error")
        print "Data, Residuals in time and assumed error"
        pylab.plot(y[0:M])
        pylab.plot(20*abs(diff[0:M]))
        pylab.plot(20*numpy.sqrt(dysq))
        
    # look at residuals in freq
    
    if display & 128:
        DIFF=Y-W*X
        pylab.figure()
        pylab.title("Data, Residuals in frequency and assumed error")    
        print "Data, Residuals in frequency and assumed error"
        pylab.plot(abs(Y))
        pylab.plot(20*abs(DIFF))
    
    #---------------------------------------------------------------------------
    # Take log and unwrap phase
    #---------------------------------------------------------------------------
    # log-amplitude and phase of W
    logW=numpy.log(W[0:N/2])
    phase_offset=0.0*numpy.pi
    "default"
    #weight=ops.good_signal(N, kstart, 0.01*N)
    "crete_data4"
    #weight=ops.good_signal(N, kstart, 0.001*N)
    weight, sstart, end = ops.good_signal_new(kstart, Y)
    
    Phasetest=ops.Phase(logW)
    logWfixedold = Phasetest.fixphase(kstart, weight, phasetolerance, phase_offset)
    logWfixed = Phasetest.newfixphase( W, kstart, weight)
    # logWfixed = Phasetest.fixphase(kstart, weight, phasetolerance, phase_offset)
    
    Wfixed = numpy.exp(logWfixed)
    Wfixedold = numpy.exp(logWfixedold)

    pylab.figure(100)
    pylab.plot(numpy.gradient(numpy.gradient(logWfixed.imag)))
    pylab.plot(weight)

    pylab.figure(2)
    pylab.xlim([-1,1])
    pylab.ylim([-1,1])

    plotRange = range(sstart,end)
    pylab.plot(W.real[plotRange], W.imag[plotRange], label = 'wrapped')
    pylab.plot(Wfixed.real[plotRange], Wfixed.imag[plotRange] , label = 'unwrapped new')
    pylab.plot(Wfixedold.real[plotRange], Wfixedold.imag[plotRange], label = 'unwrapped old')
    
#    pylab.plot(W.imag, label = 'wrapped')
#    pylab.plot(Wfixed.imag , label = 'unwrapped new')
#    pylab.plot(Wfixedold.imag, label = 'unwrapped old')

    pylab.legend()

    #pylab.show()
    
        
    
    Wnew=numpy.exp(logWfixed)
    if display & 64:
        pylab.figure()
        pylab.title("Log amplitude and phase of transfer function")
        print "Log amplitude and phase of transfer function"
        pylab.plot(logW.real)
        pylab.plot(0.01*logWfixed.imag)
        pylab.plot(0.01*logW.imag)
        pylab.plot(weight[0:N/2])
    
    #--------------------------------------------------------------------------
    # Fit model to get refractive index
    #--------------------------------------------------------------------------
    # use initial guess at thickness
    thickscale=1.0
    
    if processing[3]==0:
        print "nothing"
    #do nothing if simulating converging beam data
    
    
    if processing[3]==1:
    # simple slab
        #model=models.slab1; deriv=models.slab1dash
        (k_lo, k_hi, n, tfdiff, tfexp)=ops.convert_to_ri(kstart, thickness, models.slab1, models.slab1dash, \
                                      logWfixed, dk, numRefl)
        print "convergence failed", k_lo, k_hi
        
        "This is here to take pw data"
        #ops.savedatafilefreq(n.real, 'ipwfocppol.txt', 1, deltaf, 1)
        #ops.savedatafilefreq(n.imag, 'epwfocppol.txt', 1, deltaf, 1)
        #ops.savedatafilefreq(logWfixed.real, 'apwfocppol.txt', 1, deltaf, 1)
        #ops.savedatafilefreq(logWfixed.imag, 'ppwfocppol.txt', 1, deltaf, 1)
         
        #1/0
        
        
    if processing[3]==2:
    # GBsandwich
        #model=models.GBsand; deriv=models.GBsanddash
        (k_lo, k_hi, n)=ops.convert_to_ri(kstart, thickness, models.GBsand, models.GBsanddash, \
                                      logWfixed, dk)
        print "convergence failed", k_lo, k_hi
    
    if processing[3]==3:
    # Filled Roll
        #model=models.filled_roll; deriv=models.filled_rolldash
        (k_lo, k_hi, n)=ops.convert_to_ri(kstart, thickness, models.filled_roll, models.filled_rolldash, \
                                      logWfixed, dk)
        print "convergence failed", k_lo, k_hi
    
    
    if processing[3]==4:
    # slab in converging beam
    # first fit parallel beam
        (k_lo, k_hi, nraw, tfdiff, tfexp)=ops.convert_to_ri(kstart, thickness, models.slab1, models.slab1dash, \
                                      logWfixed, dk, numRefl)
        print "planar convergence failed", k_lo, k_hi,N/2,thickness, k_hi/paddingfactor,M/2
    
        print "raw fit"
        if display & 32 :
            pylab.figure()
            pylab.plot(nraw.real, label='index')
            pylab.plot(nraw.imag, label='extinction coefficient')
            pylab.xlabel('Frequency [THz]')
            pylab.ylabel('Amplitude [a.u.]')
            pylab.title('Raw fit')
            pylab.grid(True)
            pylab.legend()
           
    
    # set up arrays for fitting converging beam
    #    set_up_converging_fit(theta_cutoff)
    # fit converging beam
        
        "number of converging beam reflections"

        numRefl_conv = 'none'  
        #numRefl_conv = 'one' 
        #numRefl_conv = 'inf'      
    
        (k_lo, k_hi, n, tfdiffconv, tfexpconv)=ops.refine_ri(kstart, convthickness, nraw, models.convslab, \
                                         models.convslabdash, W, dk, polarisation_configuration, numRefl_conv)
        print "refinement convergence failed", k_lo, k_hi
    
        #graph(nraw[0:200].real)
        #graph(nrefined[0:200].real)
    
        #(k_lo, k_hi, nrefined)=refine_ri(kstart, thickness, nraw+0.2, convslab, \
       #                                  convslabdash, W, dk)
      #  print "refinement convergence failed", k_lo, k_hi
    
        #graph(nrefined[0:200].real)
    
    if processing[3]==5:
        "electric field parallel to the plane of incidence"
        "angle in degrees"
        angle=0.0
        sampleangle=angle*2*numpy.pi/180.0
        #print sampleangle
       # simple slab
        model=slab1para; deriv=slab1dashpara
        (k_lo, k_hi, n, tfdiff, tfexp)=ops.convert_to_ri(kstart, thickness, models.slab1para, models.slab1dashpara, \
                                      logWfixed, dk)
        print "convergence failed", k_lo, k_hi,thickness
        
    if processing[3]==6:
        "electric field perpendicular to the plane of incidence"
        "angle in degrees"
        angle=0.0
        sampleangle=angle*2*numpy.pi/180.0
        #print sampleangle
       # simple slab
        model=slab1perp; deriv=slab1dashperp
        (k_lo, k_hi, n, tfdiff, tfexp)=ops.convert_to_ri(kstart, thickness, models.slab1perp, models.slab1dashperp, \
                                      logWfixed, dk)
        print "convergence failed", k_lo, k_hi,thickness
        
    if processing[3]==7:
    # slab in converging beam
    # first fit parallel beam
        (k_lo, k_hi, nraw, tfdiff, tfexp)=ops.convert_to_ri(kstart, thickness, models.slab1, models.slab1dash, \
                                      logWfixed, dk)
        print "planar convergence failed", k_lo, k_hi
    
        print "raw fit"
        if display & 32 :
            pylab.figure()
            pylab.plot(nraw.real, label='index')
            pylab.plot(nraw.imag, label='extinction coefficient')
            pylab.xlabel('Frequency [THz]')
            pylab.ylabel('Amplitude [a.u.]')
            pylab.title('Raw fit')
            pylab.grid(True)
            pylab.legend()
           
    
    # set up arrays for fitting converging beam
    #    set_up_converging_fit(theta_cutoff)
    # fit converging beam
        (k_lo, k_hi, n, tfdiffconv, tfexpconv, fdbeamprofile)=ops.refine_ri2(kstart, convthickness, nraw, models.convslabfd, \
                                         models.convslabdashfd, W, dk)
        print "refinement convergence failed", k_lo, k_hi
        
        ops.savedatafilefreq(fdbeamprofile.real, 'beamprofile.txt', 1, deltaf)
    
        #graph(nraw[0:200].real)
        #graph(nrefined[0:200].real)
    
        #(k_lo, k_hi, nrefined)=refine_ri(kstart, thickness, nraw+0.2, convslab, \
       #                                  convslabdash, W, dk)
      #  print "refinement convergence failed", k_lo, k_hi
    
        #graph(nrefined[0:200].real)
              
    if processing[3]==8:
    # slab in converging beam donut profile
    # first fit parallel beam
        (k_lo, k_hi, nraw, tfdiff, tfexp)=ops.convert_to_ri(kstart, thickness, models.slab1, models.slab1dash, \
                                      logWfixed, dk, polarisation_configuration)
        print "planar convergence failed", k_lo, k_hi,N/2,thickness, k_hi/paddingfactor,M/2
    
        print "raw fit"
        if display & 32 :
            pylab.figure()
            pylab.plot(nraw.real, label='index')
            pylab.plot(nraw.imag, label='extinction coefficient')
            pylab.xlabel('Frequency [THz]')
            pylab.ylabel('Amplitude [a.u.]')
            pylab.title('Raw fit')
            pylab.grid(True)
            pylab.legend()
           
    
    # set up arrays for fitting converging beam
    #    set_up_converging_fit(theta_cutoff)
    # fit converging beam
        (k_lo, k_hi, n, tfdiffconv, tfexpconv)=ops.refine_ri(kstart, convthickness, nraw, models.convslabdonut, \
                                         models.convslabdashdonut, W, dk, polarisation_configuration)
        print "refinement convergence failed", k_lo, k_hi
        
    if processing[3]==9:
    # thin film on slide
        nsubindex=ops.readdatafile2("Pb_data/i_pf2_SiO2Pb00nm_R.txt", 1.0)
        nsubextin=ops.readdatafile2("Pb_data/e_pf2_SiO2Pb00nm_R.txt", 1.0)
        nnsubstrate=nsubindex+(1.0J)*nsubextin
        #ns=5.6
        Ls=1.0e-3
        (k_lo, k_hi, n, tfdiff, tfexp)=ops.convert_to_ri_thinfilm(kstart, thickness, models.thinfilm1, models.thinfilm1dash, \
                                      logWfixed, dk, nnsubstrate, Ls)
        print "planar convergence failed", k_lo, k_hi,N/2,thickness, k_hi/paddingfactor,M/2
    
        print "raw fit"
        if display & 32 :
            pylab.figure()
            pylab.plot(nraw.real, label='index')
            pylab.plot(nraw.imag, label='extinction coefficient')
            pylab.xlabel('Frequency [THz]')
            pylab.ylabel('Amplitude [a.u.]')
            pylab.title('Raw fit')
            pylab.grid(True)
            pylab.legend()  
    
    # display results using given thickness
    if display & 32:
    # magnitude of spectrum and fitted real and imag parts of refractive index
    # imag part with flipped sign to fit on graph and scaled up 10x
        print "Amplitude of transfer function, real ri, imag, ri and good window"
    # set W[0] to fix vertical scale
       
        P=4; Q=500
        pylab.figure()
        pylab.plot(fscale[P:Q], n[P:Q].real, label='conv')
        #pylab.plot(fscale[P:Q], nraw[P:Q].real, label='plane')
        pylab.title('Refractive index')
        pylab.grid(True)
        pylab.legend()
        pylab.figure()
        pylab.plot(fscale[P:Q], n[P:Q].imag, label='conv')
        #pylab.plot(fscale[P:Q], nraw[P:Q].imag, label='plane')
        pylab.title('Extinction coefficient')
        pylab.grid(True)  
        pylab.legend()  
        pylab.figure()
        pylab.plot(fscale[P:Q], abs(W[P:Q]), label='modulus of transfer function')
        pylab.plot(fscale[P:Q], weight[P:Q], label='weighting function for good signal')
        pylab.title('transfer fn + good fit')
        pylab.grid(True)
        pylab.legend()
        
    
    #----------------------------------------------------------------------------
    # Further processing
    #---------------------------------------------------------------------------
    # Fit ri to raw data
    if processing[4]==1:
    #    logWraw=log(Y/X)
    #    (k_lo, k_hi, nrefined)=refine_ri(kstart, thickness, n, model, deriv, \
    #                              logWraw, dk)
        Wraw=Y/X
        (k_lo, k_hi, nrefined)=ops.refine_ri(kstart, thickness, n, models.rawslab1, models.rawslab1dash, \
                                         Wraw, dk, polarisation_configuration)
        print "refinement convergence failed", k_lo, k_hi
 
    
    # Determine optimum thickness
    if processing[4]==2:
        lw=logW.real
        mask=ops.detrend(lw, weight, M)
        best=1e10; bestthick=0.0
        for thickscale in numpy.arange(0.9, 1.1, 0.01):
            (k_lo, k_hi, n)=ops.convert_to_ri(kstart, thickscale, models.slab1, models.slab1dash, logWfixed, kl)
            target=ops.detrend(n.imag, weight, M)
            quality=sum(weight[0:M]*target*mask)
            #print thickscale, quality
            if quality<best: best=quality; bestthick=thickscale
    
        print "best thickness", bestthick, best
        thickscale=bestthick
        (k_lo, k_hi, n)=ops.convert_to_ri(kstart, thickscale, models.slab1, models.slab1dash, logWfixed, kl)
        print "convergence failed", k_lo, k_hi
        
    # Beam profile analytical
    if processing[4]==3:
        
        "Cambridge data from 13-08-2013"
        #focusref_file = "Cambridge_data/si_at_fp1_r.txt"
        #focus_file = "Cambridge_data/si_at_fp1.txt"
        
        #focusref_file = "Cambridge_data/si_before_fp1_r.txt"
        #focus_file = "Cambridge_data/si_before_fp1.txt"
        
        #reference_file = "Cambridge_data/si_after_fp1_r.txt"
        #focus_file = "Cambridge_data/si_after_fp1.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_fully_opened1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_fully_opened1.txt"
        
        #focusref_file = "Cambridge_data/si_coll_part_first_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_first_pos1.txt"
        
        #reference_file = "Cambridge_data/nothing1.txt"
        
        "Cambridge data from 14-08-2013"
        
        #focusref_file  = "Cambridge_data/nothing2.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_second_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_second_pos1.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_third_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_third_pos1.txt"
        
        #focusref_file  = "Cambridge_data/topas_before_fp1_r.txt"
        #focus_file = "Cambridge_data/topas_before_fp1.txt"
        
        #focusref_file  = "Cambridge_data/topas_fp1_r.txt"
        #focus_file = "Cambridge_data/topas_fp1.txt"
        
        #focusref_file = "Cambridge_data/topas_coll_part1_r.txt"
        #focus_file = "Cambridge_data/topas_coll_part1.txt"
        
        #focusref_file = "Cambridge_data/3si_coll_part1_r.txt"
        #focus_file = "Cambridge_data/3si_coll_part1.txt"
        
        #focusref_file = "Cambridge_data/3si_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_fp1.txt"
        
        #focusref_file = "Cambridge_data/3si_otherside_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_otherside_fp1.txt"
        
        #focusref_file = "Cambridge_data/3si_before_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_before_fp1.txt"
        
        
        
        #focus_file = "crete_data/HDPE_focus_avg.txt"
        #focus_file = "crete_data/HDPE_collimated_00.dat"
        #focus_file = "crete_data/HDPE_focus_02.dat"
        #focus_file = "crete_data/HDPE_focus_+1mm_02.dat"
        #focus_file = "crete_data/HDPE_focus_-1mm_01.dat"
        #focus_file = "crete_data/HDPE_after_filament_02.dat"
        
        #focus_file = "Si_GaAs/gaase.txt"
        #focus_file = "Si_GaAs/gaasrec.txt"
        
      
        "Crete data2"

        #focus_file = "crete_data2/35_focus_00.dat"
        #focus_file = "crete_data2/57_focus_00.dat"  
        #focus_file = "crete_data2/35_focus_avg.txt" 
        #focus_file = "crete_data2/57_focus_avg.txt" 
        
        "Crete data4"
        #reference_file = "crete_data4/ZnTe_ref_ppol.dat"
        #sample_file = "crete_data4/HDPE_3mm_ppol_col.dat"
        focus_file = "creteSI/Si_wafer_focus_ppol_1.DAT"
        
        #reference_file = "crete_data4/ZnTe_ref_spol.dat"
        #sample_file = "crete_data4/HDPE_3mm_spol_col.dat"
        #focus_file = "crete_data4/HDPE_3mm_spol_focus.dat"
        
        
        #"mx shifts the sample scan by a number of units"
        mxf=ops.findtimeshift(reference_file, focus_file)
        #mxf=0
        
        yyf=ops.readdatafile(focus_file, 1.0, paddingfactor)
        yf=numpy.zeros(N,float)    
        yf[mxf:M+mxf]=yyf[0:M]
        #print numpy.argmax(yf)
        #yf=ops.truncate(yf, truncateval, falloff)
        Yf=fft(yf)/N; Tratio=Yf/Y; F=1.0/Tratio - 1.0; Fthickcone=1.0/numpy.sqrt(Tratio) - 1.0;
          
        theta0=numpy.zeros(N, complex)
        # gaussian
       # for k in range(2,N/2):
       #     theta0[k]=numpy.sqrt(F[k]/((0.5J)*k*dk*thickness[0]*(1.0 - 1.0/n[k])))
            
        # thin cone
        for k in range(2,N/2):
            theta0[k]=numpy.sqrt(numpy.log(1.0/Tratio[k])/((0.5J)*k*dk*thickness[0]*(1.0 - 1.0/n[k])))
            
        # thick cone
        #for k in range(2,N/2):
        #     theta0[k]=numpy.sqrt(Fthickcone[k]/((0.5J)*k*dk*thickness[0]*(1.0 - 1.0/n[k])))
            
        
        ops.savedatafilefreq(theta0.real, 'r57oatc_fpavg.txt', 1, deltaf)
        ops.savedatafilefreq(theta0.imag, 'i57oatc_fpavg.txt', 1, deltaf) 
           
        #ops.savedatafilefreq(theta0.real, 'rgasi_fp.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0.imag, 'igasi_fp.txt', 1, deltaf)

        pylab.figure()
        pylab.plot(fscale, theta0.real)
        pylab.plot(fscale, theta0.imag)
        pylab.title('theta')
        pylab.legend()
        pylab.grid(True)
        #pylab.show()
        
        # Beam profile numerical 
    if processing[4]==4:
        
        "Cambridge data from 13-08-2013"
        #focusref_file = "Cambridge_data/si_at_fp1_r.txt"
        #focus_file = "Cambridge_data/si_at_fp1.txt"
        
        #focusref_file = "Cambridge_data/si_before_fp1_r.txt"
        #focus_file = "Cambridge_data/si_before_fp1.txt"
        
        #focusref_file = "Cambridge_data/si_after_fp1_r.txt"
        #focus_file = "Cambridge_data/si_after_fp1.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_fully_opened1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_fully_opened1.txt"
        
        #focusref_file = "Cambridge_data/si_coll_part_first_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_first_pos1.txt"
        
        #reference_file = "Cambridge_data/nothing1.txt"
        
        "Cambridge data from 14-08-2013"
        
        #focusref_file  = "Cambridge_data/nothing2.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_second_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_second_pos1.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_third_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_third_pos1.txt"
        
        #focusref_file  = "Cambridge_data/topas_before_fp1_r.txt"
        #focus_file = "Cambridge_data/topas_before_fp1.txt"
        
        #focusref_file  = "Cambridge_data/topas_fp1_r.txt"
        #focus_file = "Cambridge_data/topas_fp1.txt"
        
        #focusref_file = "Cambridge_data/topas_coll_part1_r.txt"
        #focus_file = "Cambridge_data/topas_coll_part1.txt"
        
        #focusref_file = "Cambridge_data/3si_coll_part1_r.txt"
        #focus_file = "Cambridge_data/3si_coll_part1.txt"
        
        #focusref_file = "Cambridge_data/3si_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_fp1.txt"
        
        #focusref_file = "Cambridge_data/3si_otherside_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_otherside_fp1.txt"
        
        #focusref_file = "Cambridge_data/3si_before_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_before_fp1.txt"

        
        #focus_file = "crete_data/HDPE_focus_avg.txt"
        #focus_file = "crete_data/HDPE_collimated_00.dat"
        #focus_file = "crete_data/HDPE_focus_00.dat"
        #focus_file = "crete_data/HDPE_focus_+1mm_00.dat"
        #focus_file = "crete_data/HDPE_focus_-1mm_00.dat"
        #focus_file = "crete_data/HDPE_after_filament_00.dat"
        
        #focus_file = "Si_GaAs/gaase.txt"
        #focus_file = "Si_GaAs/gaasrec.txt"
        
        #focus_file = "Si/nlravg.txt"
        
        "Crete data2"

        #focus_file = "crete_data2/35_focus_00.dat"
        #focus_file = "crete_data2/57_focus_00.dat"           
        
        #focus_file = "crete_data2/35_focus_avg.txt"
        #focus_file = "crete_data2/57_focus_avg.txt"                  
                
        "Crete data3"
        
        #reference_file = "crete_data3/ZnTe_ref_ppol.dat"
        #sample_file = "crete_data3/HDPE_col_ppol.dat"
        #focus_file = "crete_data3/HDPE_focus_ppol.dat"
        
        #reference_file = "crete_data3/ZnTe_ref_spol.dat"
        #sample_file = "crete_data3/HDPE_col_spol.dat"
        #focus_file = "crete_data3/HDPE_focus_spol.dat" 
        
        "Crete data4"
        #reference_file = "crete_data4/ZnTe_ref_ppol.dat"
        #sample_file = "crete_data4/HDPE_3mm_ppol_col.dat"
        #focus_file = "crete_data4/HDPE_3mm_spol_col.dat"
        
        #reference_file = "crete_data4/ZnTe_ref_spol.dat"
        #sample_file = "crete_data4/HDPE_3mm_spol_col.dat"
        #focus_file = "samples/Teflon_converging_beam/foce.dat"
        # focus_file = "samples/HDPE/foc32.dat"
        #focus_file = "samples/oldHDPE/HDPE_focus_00.dat"

        focus_file = "samples/afoct.dat"

        pickfile ="axel.pic" 


      
        "mx shifts the sample scan by a number of units, use this for the Crete data"
        mxf=ops.findtimeshift(reference_file, focus_file)
        #mxf=0        
        
        yyf=ops.readdatafile(focus_file, 1.0, paddingfactor)
        yf=numpy.zeros(N,float)   
        
                
        yf[mxf:M+mxf]=yyf[0:M]
        
        Q=N/2
        pylab.figure()
        pylab.subplot(211)
        pylab.plot( yf[0:Q], label='yf')
        pylab.plot( x[0:Q],label='x')        
        pylab.title('scans')
        pylab.legend()
        pylab.grid(True)
        
        "crete_data4"
        #yf=ops.truncate(yf, truncateshift, falloff)
        
        
        pylab.subplot(212)
        pylab.plot( yf[0:Q], label='yf')
        pylab.plot( x[0:Q],label='x')  
        pylab.plot( y[0:Q],label='x')       
        pylab.title('truncated scans')
        pylab.legend()
        pylab.grid(True)

#        pylab.show()
        
        "default"
        #yf=ops.truncate(yf, truncateval, 30)
        
        Yf=fft(yf)/N
        
        "If a reference focus file is available then use this"
        
        #xxf=ops.readdatafile(focusref_file, 1.0, paddingfactor)
        #xf=numpy.zeros(N,float)    
        #xf[mxf:M+mxf]=xxf[0:M]
        #xf=ops.truncate(xf, truncateval, 30)
        #Xf=fft(xf)/N
        
        #theta_c=0.10
        ##theta_c=0.10+0.00J
        ##yf=models.simulate_converging_conical(reference_file, 1.0, 1.500-0.005J, 3.000e-3, theta_c, dk, 0.0, paddingfactor, N, M)       
            
        "Use reference of the parallel results"
        Wf=Yf/X
        
        "Use reference of the focus results"
        #Wf=Yf/Xf    
           
        #logWf=numpy.log(Wf)
        #phase_offset=0.0*numpy.pi
        #weightf=ops.good_signal(N, kstart, 0.01*N)
        #Phasetestf=ops.Phase(logWf)
        #logWffixed=Phasetestf.fixphase(kstart, weightf, phasetolerance, phase_offset) 
        
    
    #####analytical
        Tratio=Yf/Y; #F=1.0/Tratio - 1.0; #Fthickcone=1.0/numpy.sqrt(Tratio) - 1.0;     
        invT = Y/Yf
        Wp = Y/X             
        
        n = numpy.concatenate((n,n[::-1]))
        n = numpy.concatenate((n,[0]))
        
#        pickfile ="teflon"+str(theta)+".pic" 

        #theta = 0.06

#        Wf = Wp * (1 - 1j*numpy.pi*fscale*1e12*thickness[0]*(1-1./n)*theta**2/3e8)

#        Wf = Wp * numpy.exp(-1j*numpy.pi*fscale*1e12*convthickness[0]*(1-1./n)*theta**2/3e8)
        
#         pylab.figure()
#         pylab.plot(Tratio.real)
#         pylab.plot(invT.real)
#         pylab.plot(abs(Tratio.real)-abs(invT.real))
# 
#         #pylab.show()
# 
#         pylab.figure()
#         pylab.plot(Tratio.imag)
#         pylab.plot(invT.imag)
#         pylab.plot(abs(Tratio.imag)-abs(invT.imag))
# 
#         #pylab.show()

        
        Tratio = Wf/Wp
        theta0a=numpy.zeros(N, complex)
        theta0ai=numpy.zeros(N, complex)

        pylab.figure()
#        pylab.plot(fscale, Tratio.imag)
#        pylab.plot(fscale, numpy.arctan2(Tratio.imag,Tratio.real))

        #pylab.show()

        #theta0th=numpy.zeros(N, complex)
        
        # gaussian
       # for k in range(2,N/2):
       #     theta0[k]=numpy.sqrt(F[k]/((0.5J)*k*dk*thickness[0]*(1.0 - 1.0/n[k])))
            
        # thin cone
        F = (1-Tratio)*3e8
        
        print thickness[0]
        for k in range(2,N/2):
#            theta0a[k]=numpy.sqrt(F[k]/(1j*numpy.pi*fscale[k]*1e12*thickness[0]*(1-1./n[k])))
            theta0ai[k]=numpy.sqrt(3e8*numpy.log(1.0/Tratio[k])/(1j*numpy.pi*fscale[k]*1e12*convthickness[0]*(1-1./n[k])))
            theta0a[k]=numpy.sqrt(3e8*numpy.log(numpy.abs(Tratio[k])/Tratio[k])/(1j*numpy.pi*fscale[k]*1e12*convthickness[0]*(1-1./n[k])))

            #theta0a[k]=numpy.sqrt(3e8*(1./Tratio[k]-numpy.abs(Tratio[k]))/(1j*numpy.pi*fscale[k]*1e12*convthickness[0]*(1-1./n[k])))


            #theta0a[k]=numpy.sqrt(numpy.log(1./Tratio[k])/(1j*numpy.pi*fscale[k]*1e12*thickness[0]*(1-1./n[k])))
            
        # thick cone
        #for k in range(2,N/2):
        #     theta0th[k]=numpy.sqrt(Fthickcone[k]/((0.5J)*k*dk*thickness[0]*(1.0 - 1.0/n[k])))
            
        
        #ops.savedatafilefreq(theta0a.real, 'dhrspoltca.txt', 1, deltaf, 90/numpy.pi)
        #ops.savedatafilefreq(theta0a.imag, 'dhispoltca.txt', 1, deltaf, 90/numpy.pi) 
           



  ##########            
              
        "Gaussian and donut fitting"                    
        #(k_lo, k_hi, theta0, tfdiffconv, tfexpconv, z0indexfit, numberOfiter)=ops.convert_to_theta(kstart, convthickness, n, models.convslabtheta, \
        #                                 models.convslabdashtheta, Wf, dk, polarisation_configuration)   
        
  
        
        "hard cone fitting"
        
        (k_lo, k_hi, theta0pm, tfdiffconv, tfexpconv, numberOfiterpm)=ops.convert_to_theta_cone(kstart, convthickness, n, models.convslabthetacone, \
                                         models.convslabdashthetacone, Wf, dk, 'PM')
        
        (k_lo, k_hi, theta0pe, tfdiffconv, tfexpconv, numberOfiterpe)=ops.convert_to_theta_cone(kstart, convthickness, n, models.convslabthetacone, \
                                         models.convslabdashthetacone, Wf, dk, 'PE')

        # (k_lo, k_hi, theta0pm, tfdiffconv, tfexpconv,z0indexfit, numberOfiterpm)=ops.convert_to_theta(kstart, convthickness, n, models.convslabthetacone, \
        #                                  models.convslabdashthetacone, Wf, dk, 'PM')
        
        # (k_lo, k_hi, theta0pe, tfdiffconv, tfexpconv, z0indexfit,numberOfiterpe)=ops.convert_to_theta(kstart, convthickness, n, models.convslabthetacone, \
        #                                  models.convslabdashthetacone, Wf, dk, 'PE')

        
        #ops.savedatafilefreq(theta0pm.real, 'dhrppolhcpm.txt', 1, deltaf, 90/numpy.pi)
        #ops.savedatafilefreq(theta0pm.imag, 'dhippolhcpm.txt', 1, deltaf, 90/numpy.pi)
        #ops.savedatafilefreq(theta0pe.real,  'dhrppolhcpe.txt', 1, deltaf, 90/numpy.pi)
        #ops.savedatafilefreq(theta0pe.imag,  'dhippolhcpe.txt', 1, deltaf, 90/numpy.pi)
         
        #ops.savedatafilefreq(n.real, 'ipwcolppol.txt', 1, deltaf, 1)
        #ops.savedatafilefreq(n.imag, 'epwcolppol.txt', 1, deltaf, 1)
        #ops.savedatafilefreq(logWfixed.real, 'apwcolppol.txt', 1, deltaf, 1)
        #ops.savedatafilefreq(logWfixed.imag, 'ppwcolppol.txt', 1, deltaf, 1)
        
        
        amppm=numpy.zeros(N,complex)
        amppe=numpy.zeros(N,complex)
        ampa=numpy.zeros(N,complex)
        
        for i in range(0,N/2):
            amppm[i]=numpy.sqrt(theta0pm[i].real**2+theta0pm[i].imag**2)
            amppe[i]=numpy.sqrt(theta0pe[i].real**2+theta0pe[i].imag**2)
            ampa[i]=numpy.sqrt(theta0a[i].real**2+theta0a[i].imag**2)
            
        #ops.savedatafilefreq(amppm.real,  'dhappolhcpm.txt', 1, deltaf, 90/numpy.pi)
        #ops.savedatafilefreq(amppe.real,  'dhappolhcpe.txt', 1, deltaf, 90/numpy.pi)
        #ops.savedatafilefreq(amppe.real,  'dhaspoltca.txt', 1, deltaf, 90/numpy.pi)
        avgdiffa=0
        avgdiffpm=0
            
        import pickle        
        with open(pickfile, 'w') as f:
            pickle.dump(n, f)
            pickle.dump(theta0pe, f)
            pickle.dump(theta0a, f)
            pickle.dump(theta0ai, f)
            pickle.dump(fscale, f)
            pickle.dump(Tratio, f)

#           pickle.dump(logWfixed, f)

        
        #ops.savedatafilefreq(z0indexfit.real, 'icgpe3si_fp.txt', 1, deltaf)
        #ops.savedatafilefreq(z0indexfit.imag, 'ecgpe3si_fp.txt', 1, deltaf)
        
                   
        #ops.savedatafilefreq(theta0.real, 'rgpmr3si_bfp.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0.imag, 'igpmr3si_bfp.txt', 1, deltaf)
    
        kpos=kstart*deltaf*1e-12
           
        P=0; Q=N/2  
        pylab.figure()
        pylab.subplot(411)
        pylab.plot(fscale[P:Q], 10/numpy.pi*weight[P:Q].real, label='weighting function')
        pylab.plot(fscale[kstart], theta0pe[kstart].real, 'ro',label='kstart = '+str(kpos))        
        #pylab.plot(fscale[P:Q], logWfixed.real[P:Q],label='amplitude')
        pylab.plot(fscale[P:Q], theta0pe[P:Q].real,label='pe')
        pylab.plot(fscale[P:Q], theta0pm[P:Q].real, label='pm')       
        pylab.plot(fscale[P:Q], theta0a[P:Q].real, label='analytical') 
        pylab.ylabel('Angular spread [degrees]')
        pylab.title('beam profile real component')
        pylab.axis([0, 4, -1, 8])
        pylab.legend()
        pylab.grid(True)
        
        pylab.subplot(412)
        pylab.plot(fscale[P:Q], 10/numpy.pi*weight[P:Q].real, label='weighting function')
        pylab.plot(fscale[kstart], 90/numpy.pi*theta0pe[kstart].imag, 'ro',label='kstart = '+str(kpos))
        #pylab.plot(fscale[P:Q], logWfixed.real[P:Q],label='amplitude')
        pylab.plot(fscale[P:Q], 90/numpy.pi*theta0pe[P:Q].imag, label='pe')
        pylab.plot(fscale[P:Q], 90/numpy.pi*theta0pm[P:Q].imag, label='pm')
        pylab.plot(fscale[P:Q], 90/numpy.pi*theta0a[P:Q].imag, label='analytical')     
        pylab.ylabel('Angular spread [degrees]')
        pylab.title('beam profile imaginary component')
        pylab.axis([0, 4, -5, 5])
        pylab.legend()
        pylab.grid(True)        
                
        pylab.subplot(413)
        pylab.plot(fscale[P:Q], numberOfiterpm[P:Q], label='pm')   
        pylab.plot(fscale[P:Q], n[P:Q].real, label='col refractive index')
        pylab.plot(fscale[P:Q], n[P:Q].imag, label='col extinction coefficient')        
        pylab.ylabel('Number')
        pylab.title('Number of Newton-Raphson interations for beam profiling')
        pylab.axis([0, 4, -1, 8])
        pylab.legend()
        pylab.grid(True)
        
        pylab.subplot(414)        
        pylab.plot(fscale[kstart], logWfixed[kstart].imag, 'ro',label='kstart = '+str(kpos)) 
        pylab.plot(fscale[P:Q], logWfixed[P:Q].imag, label='unwrappedphase')
        pylab.xlabel('Frequency [THz]')
        pylab.ylabel('Phase [radians]')
        pylab.title('Unwrapped phase of collimated data')
        pylab.axis([0, 4, -120, 2])
        pylab.legend()
        pylab.grid(True)
#         
#         pylab.figure()
#         pylab.plot(fscale[P:Q], tfexpconv[P:Q].real)
#         pylab.plot(fscale[P:Q], tfexpconv[P:Q].imag)
#         pylab.title(' tfexpconv')
#         pylab.legend()
#         pylab.grid(True)
#         
#         pylab.figure()
#         pylab.plot(fscale[P:Q],  tfdiffconv[P:Q].real)
#         pylab.plot(fscale[P:Q],  tfdiffconv[P:Q].imag)
#         pylab.title('tfdiffconv')
#         pylab.legend()
#         pylab.grid(True)
#         
#         pylab.figure()
#         pylab.plot(fscale[P:Q], numberOfiter[P:Q])
#         pylab.title('numberOfiter')
#         pylab.legend()
#         pylab.grid(True)
#         
#         pylab.figure()
#         pylab.plot(fscale[P:Q], numpy.log(tfexpconv[P:Q]).real, label='extracted')
#         pylab.plot(fscale[P:Q], numpy.log(Wf[P:Q]).real, label='experimental')
#         pylab.title('tfconv amp')
#         pylab.legend()
#         pylab.grid(True)
#         
#         pylab.figure()
#         pylab.plot(fscale[P:Q], numpy.log(tfexpconv[P:Q]).imag, label='extracted')
#         pylab.plot(fscale[P:Q], numpy.log(Wf[P:Q]).imag, label='experimental')
#         pylab.title('tfconv phase')
#         pylab.legend()
#         pylab.grid(True)
        
        pylab.show()
        
    if processing[4]==5:
        
        polarisation_configuration_simulated = 'PE'
        numRefl_simulated = 'none'
        #numRefl_simulated = 'inf'        
        
        numRefl_simulated_plane = 'none'   
        #numRefl_simulated_plane = 'inf'             
        layers=[(3.500-0.000J, 0.525e-3)]
        
        ya=models.simulate_parallel2(reference_file, layers, dk, 0.0, paddingfactor, N, M, numRefl_simulated_plane)
        Ya=fft(ya)/N;
           
        angular_spread=numpy.zeros(N/2, complex)
        for k in range(0,N/2):
            #angular_spread[k]= 0.7/(1.0 + numpy.exp(float(k)/100))+0.2/(1.0 + float(k)/N)
            #angular_spread[k]= 0.7/(1.0 + numpy.exp(float(k)/100))+0.2/(1.0 + float(k)/N)
            angular_spread[k]=angular_spread_sim
        
        P=0; Q=N/8  
        pylab.figure()
        pylab.plot(fscale[P:Q], angular_spread[P:Q].real, label="real") 
        pylab.plot(fscale[P:Q], angular_spread[P:Q].imag, label="imag") 
        pylab.legend()
        #pylab.show()
        #yf=models.simulate_converging_conical(reference_file, 1.0, 1.500-0.000J, 3.000e-3, 0.1, dk, 0.0, paddingfactor, N, M) 
        
        yf=models.simulate_converging(reference_file, 1.0, 3.500-0.000J, 0.525e-3, angular_spread, dk, 0.0, paddingfactor, N, M, polarisation_configuration_simulated, numRefl_simulated)      
        
        Yf=fft(yf)/N
        Wf=Yf/X    
        nsim=numpy.zeros(N,complex)
        nsim[0:N]=3.500-0.000J        
        
        "Gaussian"
                                                   
        (k_lo, k_hi, theta0pm, tfdiffconv, tfexpconv, z0indexfit, numberOfiter)=ops.convert_to_theta(kstart, convthickness, nsim, models.convslabtheta, \
                                         models.convslabdashtheta, Wf, dk, 'PM')        
                                                           
        (k_lo, k_hi, theta0pe, tfdiffconv, tfexpconv, z0indexfit, numberOfiter)=ops.convert_to_theta(kstart, convthickness, nsim, models.convslabtheta, \
                                         models.convslabdashtheta, Wf, dk, 'PE')
        
        #hard cone
        
        #(k_lo, k_hi, theta0, tfdiffconv, tfexpconv)=ops.convert_to_theta(kstart, convthickness, nsim, models.convslabthetacone, \
        #                                 models.convslabdashthetacone, Wf, dk, polarisation_configuration)
        
        
        

        
        
        #analytical method
     
        Yfa=fft(yf)/N; Tratio=Yfa/Ya
        
        F=1.0/Tratio - 1.0 
        theta0a=numpy.zeros(N, complex)
        # gaussian
        
        for k in range(0,N/2):
            theta0a[k]=numpy.sqrt(F[k]/((0.5J)*k*dk*thickness[0]*(1.0 - 1.0/nsim[k])))  
            
        #difference
        diffa = numpy.zeros(N, complex)
        diffpm = numpy.zeros(N, complex)
        for k in range(0,N/2):
            diffa[k] = theta0pe[k] - theta0a[k] 
            diffpm[k] = theta0pe[k] - theta0pm[k] 
            
        
        #k=435 = 1 THz, k = 1738 = 4 THz for 8x padding of Cambridge data
        
        avgdiffa = numpy.mean(diffa[435:1738]) 
        avgdiffpm = numpy.mean(diffpm[435:1738])  
            
        #ops.savedatafilefreq(theta0a.real, 'rsgpega.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0a.imag, 'isgpega.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0pm.real, 'rsgpegpm.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0pm.imag, 'isgpegpm.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0pe.real, 'rsgpegpe.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0pe.imag, 'isgpegpe.txt', 1, deltaf)
        
        #ops.savedatafilefreq(diffa.real, 'darsgpe.txt', 1, deltaf)
        #ops.savedatafilefreq(diffa.imag, 'daisgpe.txt', 1, deltaf)
        #ops.savedatafilefreq(diffpm.real, 'dpmrsgpe.txt', 1, deltaf)
        #ops.savedatafilefreq(diffpm.imag, 'dpmisgpe.txt', 1, deltaf)
        
                
        P=0; Q=N/2   
       # pylab.figure()
       # pylab.plot(fscale[P:Q], tfdiffconv[P:Q].real, label="N")
       # pylab.plot(fscale[P:Q], tfdiffconv[P:Q].imag, label="N")
       # pylab.title('difference')
       # pylab.legend()
       # pylab.grid(True)
        
       # pylab.figure()
       # pylab.plot(fscale[P:Q], tfexpconv[P:Q].real)
       # pylab.plot(fscale[P:Q], Wf[P:Q].real)
       # pylab.title('tfamplitude')
       # pylab.legend()
       # pylab.grid(True)
        
       # pylab.figure()
       # pylab.plot(fscale[P:Q], tfexpconv[P:Q].imag)
       # pylab.plot(fscale[P:Q], Wf[P:Q].imag)
        #pylab.title('tfphase')
        #pylab.legend()
        #pylab.grid(True)
                
        P=0; Q=N/8   
        pylab.figure()
        pylab.plot(fscale[P:Q], diffa[P:Q].real, label="dar")
        pylab.plot(fscale[P:Q], diffa[P:Q].imag, label="dai")
        pylab.plot(fscale[P:Q], diffpm[P:Q].real, label="dpmr")
        pylab.plot(fscale[P:Q], diffpm[P:Q].imag, label="dpmi")
        pylab.legend()
        
        pylab.figure()
        pylab.plot(fscale[P:Q], theta0a[P:Q].real, label="ar")
        pylab.plot(fscale[P:Q], theta0a[P:Q].imag, label="ai")
        pylab.plot(fscale[P:Q], theta0pe[P:Q].real, label="per")
        pylab.plot(fscale[P:Q], theta0pe[P:Q].imag, label="pei")
        pylab.plot(fscale[P:Q], theta0pm[P:Q].real, label="pmr")
        pylab.plot(fscale[P:Q], theta0pm[P:Q].imag, label="pmi")
        pylab.title('theta')
        pylab.legend()
        pylab.grid(True)
        #pylab.show()
        
    if processing[4]==6:
        
        "Cambridge data from 13-08-2013"
        #focusref_file = "Cambridge_data/si_at_fp1_r.txt"
        #focus_file = "Cambridge_data/si_at_fp1.txt"
        
        #focusref_file = "Cambridge_data/si_before_fp1_r.txt"
        #focus_file = "Cambridge_data/si_before_fp1.txt"
        
        #focusref_file = "Cambridge_data/si_after_fp1_r.txt"
        #focus_file = "Cambridge_data/si_after_fp1.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_fully_opened1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_fully_opened1.txt"
        
        #focusref_file = "Cambridge_data/si_coll_part_first_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_first_pos1.txt"
        
        #reference_file = "Cambridge_data/nothing1.txt"
        
        "Cambridge data from 14-08-2013"
        
        #focusref_file  = "Cambridge_data/nothing2.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_second_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_second_pos1.txt"
        
        #focusref_file  = "Cambridge_data/si_coll_part_third_pos1_r.txt"
        #focus_file = "Cambridge_data/si_coll_part_third_pos1.txt"
        
        #focusref_file  = "Cambridge_data/topas_before_fp1_r.txt"
        #focus_file = "Cambridge_data/topas_before_fp1.txt"
        
        #focusref_file  = "Cambridge_data/topas_fp1_r.txt"
        #focus_file = "Cambridge_data/topas_fp1.txt"
        
        #focusref_file = "Cambridge_data/topas_coll_part1_r.txt"
        #focus_file = "Cambridge_data/topas_coll_part1.txt"
        
        #focusref_file = "Cambridge_data/3si_coll_part1_r.txt"
        #focus_file = "Cambridge_data/3si_coll_part1.txt"
        
        #focusref_file = "Cambridge_data/3si_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_fp1.txt"
        
        #focusref_file = "Cambridge_data/3si_otherside_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_otherside_fp1.txt"
        
        #focusref_file = "Cambridge_data/3si_before_fp1_r.txt"
        #focus_file = "Cambridge_data/3si_before_fp1.txt"

        
        #focus_file = "crete_data/HDPE_collimated_00.dat"
        #focus_file = "crete_data/HDPE_focus_00.dat"
        #focus_file = "crete_data/HDPE_focus_+1mm_00.dat"
        #focus_file = "crete_data/HDPE_focus_-1mm_00.dat"
        #focus_file = "crete_data/HDPE_after_filament_00.dat"
        
        #focus_file = "Si_GaAs/gaase.txt"
        #focus_file = "Si_GaAs/gaasrec.txt"
        
        #focus_file = "Si/nlravg.txt"
                
        "mx shifts the sample scan by a number of units"
        #mxf=ops.findtimeshift(reference_file, focus_file)
        mxf=0
        
        
        yyf=ops.readdatafile(focus_file, 1.0, paddingfactor)
        yf=numpy.zeros(N,float)    
        yf[mxf:M+mxf]=yyf[0:M]
        #yf=ops.truncate(yf, truncateval, 30)
        Yf=fft(yf)/N
        
        "If a reference focus file is available then use this"
        
        #xxf=ops.readdatafile(focusref_file, 1.0, paddingfactor)
        #xf=numpy.zeros(N,float)    
        #xf[mxf:M+mxf]=xxf[0:M]
        #xf=ops.truncate(xf, truncateval, 30)
        #Xf=fft(xf)/N
        
        #theta_c=0.10
        ##theta_c=0.10+0.00J
        ##yf=models.simulate_converging_conical(reference_file, 1.0, 1.500-0.005J, 3.000e-3, theta_c, dk, 0.0, paddingfactor, N, M)       
            
        "Use reference of the parallel results"
        Wf=Yf/X
        
        "Use reference of the focus results"
         #Wf=Yf/Xf            
                
              
        "Gaussian and donut fitting with least sqaures fitting"                    
        (k_lo, k_hi, theta0, tfdiffconv, tfexpconv)=ops.convert_to_theta_leastsq(kstart, convthickness, n, models.convslabtheta, \
                                         models.convslabdashtheta, Wf, dk, polarisation_configuration)         
                   
        #ops.savedatafilefreq(theta0.real, 'rgpmrsi_c2.txt', 1, deltaf)
        #ops.savedatafilefreq(theta0.imag, 'igpmrsi_c2.txt', 1, deltaf)
        
        P=0; Q=N/2   
        pylab.figure()
        pylab.plot(fscale[P:Q], theta0[P:Q].real)
        #pylab.plot(fscale[P:Q], theta0[P:Q].imag)
        pylab.title('theta')
        pylab.legend()
        pylab.grid(True)
        #pylab.show()

    if processing[4]==7:
        
        #Step 1. Determine the bp from the Cambridge data
        
        "Cambridge data from 14-08-2013"        
        focus_file = "afoc.dat"   
        
        'Create converging beam transfer function using reference from parallel results'     
        mxf=0       
        yyf=ops.readdatafile(focus_file, 1.0, paddingfactor)
        yf=numpy.zeros(N,float)    
        yf[mxf:M+mxf]=yyf[0:M]
        yfold = yf[:]
        pylab.figure(10)
        pylab.plot(yf)
#        yf=ops.truncate(yf, truncateshift, falloff)
        pylab.plot(yf)
        pylab.plot(yf/yfold)
        #pylab.show()
        Yf=fft(yf)/N            
        Wf=Yf/X        
        
        (k_lo, k_hi, theta0realpm, tfdiffconv, tfexpconv, z0indexfit, numberOfiter)=ops.convert_to_theta(kstart, convthickness, n, models.convslabtheta, \
                                         models.convslabdashtheta, Wf, dk, polarisation_configuration)   
          
        #Step2. Simulate the data using the bp
        
        polarisation_configuration_simulated = 'PE'
        numRefl_simulated = 'none'
        #numRefl_simulated = 'inf'        
        
        numRefl_simulated_plane = 'none'   
        #numRefl_simulated_plane = 'inf'         
        
        yfsim=models.simulate_converging_freq(reference_file, 1.0, n, convthickness, theta0realpm, dk, 0.0, paddingfactor, N, M, polarisation_configuration_simulated, numRefl_simulated)      
               
        Yfsim=fft(yfsim)/N
        Wfsim=Yfsim/X      
        
        print "simulated data using beam profile"   
        
        #Step3. Determine the bp using the different methods 
        
        "Gaussian"
                                                   
        (k_lo, k_hi, theta0pm, tfdiffconv, tfexpconv, z0indexfit, numberOfiter)=ops.convert_to_theta(kstart, convthickness, n, models.convslabtheta, \
                                        models.convslabdashtheta, Wfsim, dk, 'PM')  
        
        print "determined numerical PM"       
                                                           
        (k_lo, k_hi, theta0pe, tfdiffconv, tfexpconv, z0indexfit, numberOfiter)=ops.convert_to_theta(kstart, convthickness, n, models.convslabtheta, \
                                         models.convslabdashtheta, Wfsim, dk, 'PE')    
         
        print "determined numerical PE"    
         
        #analytical method
      
        Tratio=Yfsim/Y; F=1.0/Tratio - 1.0 
        theta0a=numpy.zeros(N, complex)
        # gaussian
         
        for k in range(0,N/2):
            theta0a[k]=numpy.sqrt(F[k]/((0.5J)*k*dk*thickness[0]*(1.0 - 1.0/n[k])))  
             
        print "determined analytical"  
             
        #difference
        diffa = numpy.zeros(N, complex)
        diffpm = numpy.zeros(N, complex)
        for k in range(0,N/2):
            diffa[k] = theta0pe[k] - theta0a[k] 
            diffpm[k] = theta0pe[k] - theta0pm[k] 
            
        
        #k=435 = 1 THz, k = 1738 = 4 THz for 8x padding of Cambridge data
        
        avgdiffa = 0
        avgdiffpm = 0
        
        #avgdiffa = numpy.mean(diffa[435:1738]) 
        #avgdiffpm = numpy.mean(diffpm[435:1738])  
            
#         ops.savedatafilefreq(theta0a.real, 'rsgpega.txt', 1, deltaf)
#         ops.savedatafilefreq(theta0a.imag, 'isgpega.txt', 1, deltaf)
#         ops.savedatafilefreq(theta0pm.real, 'rsgpegpm.txt', 1, deltaf)
#         ops.savedatafilefreq(theta0pm.imag, 'isgpegpm.txt', 1, deltaf)
#         ops.savedatafilefreq(theta0pe.real, 'rsgpegpe.txt', 1, deltaf)
#         ops.savedatafilefreq(theta0pe.imag, 'isgpegpe.txt', 1, deltaf)
#         
#         ops.savedatafilefreq(diffa.real, 'darsgpe.txt', 1, deltaf)
#         ops.savedatafilefreq(diffa.imag, 'daisgpe.txt', 1, deltaf)
#         ops.savedatafilefreq(diffpm.real, 'dpmrsgpe.txt', 1, deltaf)
#         ops.savedatafilefreq(diffpm.imag, 'dpmisgpe.txt', 1, deltaf)
#         
                
        P=0; Q=N/2   
       # pylab.figure()
       # pylab.plot(fscale[P:Q], tfdiffconv[P:Q].real, label="N")
       # pylab.plot(fscale[P:Q], tfdiffconv[P:Q].imag, label="N")
       # pylab.title('difference')
       # pylab.legend()
       # pylab.grid(True)
        
       # pylab.figure()
       # pylab.plot(fscale[P:Q], tfexpconv[P:Q].real)
       # pylab.plot(fscale[P:Q], Wf[P:Q].real)
       # pylab.title('tfamplitude')
       # pylab.legend()
       # pylab.grid(True)
        
       # pylab.figure()
       # pylab.plot(fscale[P:Q], tfexpconv[P:Q].imag)
       # pylab.plot(fscale[P:Q], Wf[P:Q].imag)
        #pylab.title('tfphase')
        #pylab.legend()
        #pylab.grid(True)
                
        import pickle
        
        with open("a"+str(truncateshift)+"_"+str(falloff)+".pic", 'w') as f:
            pickle.dump(n, f)
            pickle.dump(theta0pe, f)
            pickle.dump(theta0a, f)
        
        
        P=0; Q=N/8   
        pylab.figure()
        pylab.plot(fscale[P:Q],n[P:Q].real, label='n real')
        pylab.plot(fscale[P:Q],n[P:Q].imag, label='n imag')
        
        pylab.figure()
        pylab.plot(fscale[P:Q], diffa[P:Q].real, label="dar")
        pylab.plot(fscale[P:Q], diffa[P:Q].imag, label="dai")
        pylab.plot(fscale[P:Q], diffpm[P:Q].real, label="dpmr")
        pylab.plot(fscale[P:Q], diffpm[P:Q].imag, label="dpmi")
        pylab.legend()
         
        pylab.figure()
        pylab.plot(fscale[P:Q], theta0a[P:Q].real, label="ar")
        pylab.plot(fscale[P:Q], theta0a[P:Q].imag, label="ai")
        pylab.plot(fscale[P:Q], theta0pe[P:Q].real, label="per")
        pylab.plot(fscale[P:Q], theta0pe[P:Q].imag, label="pei")
        pylab.plot(fscale[P:Q], theta0pm[P:Q].real, label="pmr")
        pylab.plot(fscale[P:Q], theta0pm[P:Q].imag, label="pmi")
        pylab.title('theta')
        pylab.legend()
        pylab.grid(True)
        #pylab.show()
        

    
    #----------------------------------------------------------------------------
    # Show results of further processing
    #----------------------------------------------------------------------------
    
    if display & 16:
        print "Futher Results"
    # set W[0] to fix vertical scale
        W[0]=3.0
    # magnitude of spectrum and fitted real and imag parts of refractive index
    # imag part with flipped sign to fit on graph and scaled up 10x
        P=0; Q=N/2
        pylab.figure()
        pylab.plot(fscale[P:Q], n[P:Q].real, label='index')
        pylab.xlabel('Frequency [THz]')
        pylab.ylabel('Amplitude [a.u.]')
        pylab.title('Further results')
        pylab.grid(True)
        pylab.plot(fscale[P:Q], n[P:Q].imag, label='extinction coefficient')
        pylab.plot(fscale[P:Q], abs(W[P:Q]), label='modulus of transfer function')
        pylab.plot(fscale[P:Q], weight[0:Q], label='weighting function for good signal')
        pylab.legend()
   
    avgdiffa=0
    avgdiffpm=0
    return (W, logW, logWfixed, n, tfdiff, tfexp, fscale, tscale, weight, Wnew, w, weight, avgdiffa, avgdiffpm)
    #return (W, logW, logWfixed, n, tfdiff, tfexp, fscale, tscale, weight, Wnew, w, weight)




































