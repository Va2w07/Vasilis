import numpy
import process
import pylab
import ops

# #fileobj=open("beamprofile.txt", 'r')
# 
# bp=[]; f=[]
# for line in fileobj.readlines():
#     l=line.split(",")
#     bp.append(complex(l[1]))
#     f.append(float(l[0]))
# fileobj.close()
# bp=numpy.array(bp)
# f=numpy.array(f)
# 
# M=len(bp)
# fileobj=open('beam.txt', 'w')
# for i in range(0,M):
#     fileobj.write(str(f[i])+","+str(bp[i].real)+"\n")
# fileobj.close() 
        

#import models
##########################################################################
# Set the main parameters of the data set
##########################################################################

    
 #---------------------------------------------------------------------------
 # Read or create data files
 #---------------------------------------------------------------------------
 
"sample and reference file names"

"Cambridge data from 13-08-2013"

#reference_file = "Cambridge_data/si_at_fp1_r.txt"
#sample_file = "Cambridge_data/si_at_fp1.txt"

#reference_file = "Cambridge_data/si_before_fp1_r.txt"
#sample_file = "Cambridge_data/si_before_fp1.txt"

#reference_file = "Cambridge_data/si_after_fp1_r.txt"
#sample_file = "Cambridge_data/si_after_fp1.txt"

#reference_file = "Cambridge_data/si_coll_part_fully_opened1_r.txt"
#sample_file = "Cambridge_data/si_coll_part_fully_opened1.txt"

#reference_file = "Cambridge_data/si_coll_part_first_pos1_r.txt"
#sample_file = "Cambridge_data/si_coll_part_first_pos1.txt"

#reference_file = "Cambridge_data/nothing1.txt"

"Cambridge data from 14-08-2013"

#reference_file = "Cambridge_data/nothing2.txt"

#reference_file = "Cambridge_data/si_coll_part_second_pos1_r.txt"
#sample_file = "Cambridge_data/si_coll_part_second_pos1.txt"

#reference_file = "Cambridge_data/si_coll_part_third_pos1_r.txt"
#sample_file = "Cambridge_data/si_coll_part_third_pos1.txt"

#reference_file = "Cambridge_data/topas_before_fp1_r.txt"
#sample_file = "Cambridge_data/topas_before_fp1.txt"

#reference_file = "Cambridge_data/topas_fp1_r.txt"
#sample_file = "Cambridge_data/topas_fp1.txt"

#reference_file = "Cambridge_data/topas_coll_part1_r.txt"
#sample_file = "Cambridge_data/topas_coll_part1.txt"

#reference_file = "Cambridge_data/3si_coll_part1_r.txt"
#sample_file = "Cambridge_data/3si_coll_part1.txt"

#reference_file = "Cambridge_data/3si_fp1_r.txt"
#sample_file = "Cambridge_data/3si_fp1.txt"

#reference_file = "Cambridge_data/3si_otherside_fp1_r.txt"
#sample_file = "Cambridge_data/3si_otherside_fp1.txt"

#reference_file = "Cambridge_data/3si_before_fp1_r.txt"
#sample_file = "Cambridge_data/3si_before_fp1.txt"



"Other Samples"
#reference_file = "Si/nlrefavg.txt"
#sample_file = "Si/nlcavg.txt"

#reference_file = "topas/topas_reference.txt"
#sample_file = "topas/topas_sample.txt"

#reference_file = "spectrosil0.5mm2.txt"
#sample_file = "spectrosil0.5mm.txt"

"Crete data"

#reference_file = "crete_data/Reference_avg.txt"
#sample_file = "crete_data/HDPE_collimated_avg.txt"
#sample_file = "crete_data/HDPE_focus_avg.txt"

#reference_file = "crete_data/Reference_00.dat"
#sample_file = "crete_data/HDPE_collimated_00.dat"
#sample_file = "crete_data/HDPE_after_filament_00.dat"
#sample_file = "crete_data/HDPE_focus_00.dat"
#sample_file = "crete_data/HDPE_focus_+1mm_00.dat"
#sample_file = "crete_data/HDPE_focus_-1mm_00.dat"

"Crete data2"
#reference_file = "crete_data2/35_Reference_00.dat"
#sample_file = "crete_data2/35_collimated_00.dat"
#sample_file = "crete_data2/35_focus_00.dat"

#reference_file = "crete_data2/35_Reference_avg.txt"
#sample_file = "crete_data2/35_collimated_avg.txt"
#sample_file = "crete_data2/35_focus_avg.txt"

#reference_file = "crete_data2/57_Reference_avg.txt"
#sample_file = "crete_data2/57_collimated_avg.txt"
#sample_file = "crete_data2/57_focus_avg.txt"

#reference_file = "crete_data2/57_Reference_00.txt"
#sample_file = "crete_data2/57_collimated_00.txt"
#sample_file = "crete_data2/57_focus_00.txt"

"Crete data3"
#reference_file = "crete_data3/ZnTe_ref_ppol.dat"
#sample_file = "crete_data3/HDPE_col_ppol.dat"
#sample_file = "crete_data3/HDPE_focus_ppol.dat"

#reference_file = "crete_data4/ZnTe_ref_spol.dat"
#sample_file = "crete_data4/HDPE_col_spol.dat"
#sample_file = "crete_data3/HDPE_focus_spol.dat"

"teflon"
# reference_file = "samples/Teflon_converging_beam/teflon.20mm.ref"
# sample_file = "samples/Teflon_converging_beam/teflon.20mm.collimated.sam"

#reference_file = "samples/Teflon_converging_beam/refe.dat"
#sample_file = "samples/Teflon_converging_beam/cole.dat"


"Crete data4"

#reference_file = "samples/HDPE/ref.dat"
#sample_file = "samples/HDPE/col.dat"

# reference_file = "samples/HDPE/ref3.dat"
# sample_file = "samples/HDPE/col3.dat"

#sample_file = "crete_data4/HDPE_3mm_spol_focus.dat"

reference_file = "samples/aref.dat"
sample_file = "samples/acolt.dat"
#sample_file = "crete_data4/HDPE_3mm_spol_focus.dat"

"PB thin film (Marks)"
#reference_file = "Pb_data/Reference.txt"
#reference_file = "Pb_data/SiO2Pb00nm.txt"
#sample_file = "Pb_data/SiO2Pb80nm.txt"
#sample_file = "Pb_data/SiO2Pb60nm.txt"
#sample_file = "Pb_data/SiO2Pb40nm.txt"
#sample_file = "Pb_data/PbSiO280nm.txt"
#sample_file = "Pb_data/SiO2Pb00nm.txt"


#reference_file = "Si_GaAs/ref.txt"
#sample_file = "Si_GaAs/gaaspara.txt"
#sample_file = "Si_GaAs/gaase3.txt"
#sample_file = "Si_GaAs/gaasrec.txt"

"thickness of sample for the planewave extraction"
#thickness = [0.477e-3] #"si"
#thickness = [2.670e-3] #"topas"
#thickness = [0.460e-3] #"3si"
thickness = [0.525e-3] #"3si"

#thickness = [21.0e-3] #"PTFE"
#thickness = [3.000e-3] #"HDPE"
# thickness = [5.000e-3] #"HDPE"

"thickness of sample for the converging beam extraction"
#convthickness = [0.477e-3] #"si"
#convthickness = [2.670e-3] #"topas"
#convthickness = [0.460e-3] #"3si"
convthickness = [0.525e-3] #"3si"
#convthickness = [21.0e-3] #"PTFE"   
#convthickness = [3.000e-3] #"HDPE"   
# convthickness = [5.000e-3] #"HDPE"   
    

"sets the converging beam extraction angle"
setkmin=20
setthetacutoff=0.3
setangularcentre=0.2

# size of data record
#M = 8096
"Cambridge data has 1999 data points"
M=1999
    
"spectrosil data has 2048 data points"
#M=2048
"crete data, creta data 2"
#M=1202

"crete data 3"
#M=451

"crete data 4"
#M=401

"Pb/Si, GaAs data, topas, Si (3mm)"
# M=702

"padding factor sets the number of multiples which the data is padded with zero, 2 doubles the data length with zeros, 3 triples the data length with zeros."
paddingfactor=5
"N is double the number of data points, in the read data file definition the data is automatically padded with zeros to double the number of points there is a m somewhere instead of N"
N=paddingfactor*M
# sampling interval 
"deltat is the time step between data points"

"Cambridge data has 1999 data points"
deltat=2.716e-14

"spectrosil data has 2048 data points"
#deltat=1.93e-14

"crete data, crete data 2, crete data 3"
#deltat=2.00e-14

"crete data 4"
# deltat=19.9412e-15

#deltat = 0.0488e-12
#deltat = 0.0053e-12
"Pb/Si, GaAs data, Si"
#deltat=13.333e-14
"topas 2.67mm"
#deltat=14.6484375e-14

"deltaf is the frequency step between data points in the frequency domain"
deltaf=1.0/(N*deltat)
# wave number for frequency point k is k*dk
"dk is the wave number step in frequency, dk allows a code efficient way to step up and down in frequency using integer steps"
dk=2*numpy.pi*deltaf/2.9979e8 

"number of plane wave reflections"

numRefl = 'none'  
#numRefl = 'one' 
#numRefl = 'inf'   

"PM or PE configuration for the converging beam"
 
#polarisation_configuration = 'PM'
polarisation_configuration = 'PE'
 
##########################################################################
# MAIN PROGRAM
############################################################################

#
# Choose what to do
# Data source:
#     sample file                   0
#     simulated parallel beam       1
#     simulated sandwich            2
#     simulated converging beam     3
#     simulated converging beam     
#     using conical beam profile    4

# Truncation:
#     none                          0
#     default setting               1

# Noise Model and Treatment:
#     constant noise no supression  0
#     constant noise supression     1
#     transparent default           2
#     constant noise with causality 3 (uses Levinson)
#     signal dependent noise        4

# Interpretation Model:
#     none if simulating converging 
#     beam data                     0
#     parallel single slab          1
#     parallel GB sandwich          2
#     parallel filled roll          3 (fits outer layers)
#     converging beam               4
#     planewave parallel angle      5
#     planewave perpendicular angle 6
#     beam profiling                7
#     donut profile                 8
#     Pb thin film                  9

# Further processing
#     none                          0
#     refit to raw data             1
#     optimise thickness            2
#     beam profile analytical       3  
#     beam profile numerical        4
 
#     beam profile      
#     simulated data                5 
#     beam profile numerical
#     with least squares fitting    6 
#     extract bp from real data then
#     simulate that data with bp and
#     extract with different 
#     methods                       7  

processing=(0, 1, 0, 1, 4)

"default"
kstart=30*paddingfactor
 
"crete_data4"   
#kstart=10*paddingfactor


phasetolerance=0.3

    
# Choose what to see.
# raw data                                               2048
# truncate the data and display                          1024
# impulse response time domain                           512
# residuals in time                                      256
# residuals in frequency                                 128
# show unwrapped phase and log amplitude                 64
# show raw ri estimates                                  32
# show refined ri estimates                              16
# everything=4095

display=0



# def phasecheck(phase):
#     g=len(phase)
#     phasediff=numpy.zeros(g)
#     for i in numpy.arange(0,g-1):
#         if i-1<0:
#             phasediff[i]=phase[i]-(phase[i+1]+phase[i])/2
#         else:
#             phasediff[i]=phase[i]-(phase[i+1]+phase[i-1])/2      
#     return phasediff
#     
# 
# phasecalc=numpy.zeros((13,N/2), complex)
# RI=numpy.zeros((13,N/2), complex)
# phasedifference=numpy.zeros((13,N/2), complex)
# weight1=numpy.zeros(M,complex)
# 
# k=0
# for i in numpy.arange(0.4,1.7,0.1):
#     (W, logW, logWfixed, n, tfdiff, tfexp, fscale, tscale, weight, Wnew, w, weight)=process.extraction(N, M, dk, deltat, deltaf, reference_file, sample_file, thickness, convthickness, \
#                    setthetacutoff, processing, kstart, setkmin, paddingfactor, i, setangularcentre, display)
#     phasecalc[k]=logWfixed.imag
#     RI[k]=n.real
#     phasedifference[k]=phasecheck(logWfixed.imag)
#     weight1=weight
#     k+=1
#     
#     
# #pylab.figure()
# #pylab.plot(phasedifference[1])
# #pylab.title('Refractive index')
# #pylab.legend()
# 
# 
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib
# from matplotlib import cm
# ax = pylab.figure().add_subplot(111, projection='3d')
# 
# 
# 
# # create supporting points in polar coordinates
# x = numpy.linspace(0.4,1.6,13)
# y = fscale[0:N/2]
# X,Y = numpy.meshgrid(x,y)
# 
# 
# #Z = phasecalc
# Z = phasedifference
# 
# indexgrad=numpy.zeros((N/2,13), complex)
# for i in numpy.arange(0,N/2):
#     indexgrad[i]=numpy.gradient(RI.T[i])   
#   
# 
# 
# ax.plot_wireframe(X, Y, Z.T, rstride=30, cstride=1)
# #ax.plot_surface(X, Y, Z.T, rstride=30, cstride=1, cmap=cm.coolwarm,
# #        linewidth=0, antialiased=False)
# #ax.set_zlim3d(0, 1)
# #ax.set_ylim3d(0, 2)
# ax.set_xlabel('phase tolerance [radians]')
# ax.set_ylabel('frequency [THz]')
# ax.set_zlabel('phase [radians]')
# #pylab.show()
# 
# #pylab.figure()
# #pylab.imshow(phasedifference.imag)
# 
# 
# bx = pylab.figure().add_subplot(111, projection='3d')
# bx.plot_wireframe(X, Y, RI.T, rstride=30, cstride=1)
# #bx.plot_surface(X, Y, RI.T, rstride=30, cstride=1, cmap=cm.coolwarm,
# #        linewidth=0, antialiased=False)
# bx.set_zlim3d(1.40 , 1.53)
# #ax.set_ylim3d(0, 2)
# bx.set_xlabel('phase tolerance [radians]')
# bx.set_ylabel('frequency [THz]')
# bx.set_zlabel('Refractive index [a.u.]')
# 
# cx = pylab.figure().add_subplot(111, projection='3d')
# cx.plot_wireframe(X, Y, indexgrad, rstride=30, cstride=1)
# #bx.plot_surface(X, Y, RI.T, rstride=30, cstride=1, cmap=cm.coolwarm,
# #        linewidth=0, antialiased=False)
# cx.set_zlim3d(-0.1 , 0.1)
# #ax.set_ylim3d(0, 2)
# cx.set_xlabel('phase tolerance [radians]')
# cx.set_ylabel('frequency [THz]')
# cx.set_zlabel('gradient between index [a.u.]')
#  
# pylab.show() 
#  
# 1/0   
# pylab.figure()    
# pylab.xlabel('Input')
# pylab.ylabel('Function values')
# pylab.plot_surface(phasecalc.imag)
# pylab.show()

angular_spread_sim = 0.0

#for theta in numpy.arange(0.0,1,0.01):
#    print theta
#    (W, logW, logWfixed, n, tfdiff, tfexp, fscale, tscale, weight, Wnew, w, weight, avgdiffa, avgdiffpm)=process.extraction(N, M, dk, deltat, deltaf, reference_file, sample_file, thickness, convthickness, \
#                   setthetacutoff, processing, kstart, setkmin, paddingfactor, phasetolerance, setangularcentre, display, numRefl, polarisation_configuration, angular_spread_sim, theta)
theta = 0
(W, logW, logWfixed, n, tfdiff, tfexp, fscale, tscale, weight, Wnew, w, weight, avgdiffa, avgdiffpm)=process.extraction(N, M, dk, deltat, deltaf, reference_file, sample_file, thickness, convthickness, \
		setthetacutoff, processing, kstart, setkmin, paddingfactor, phasetolerance, setangularcentre, display, numRefl, polarisation_configuration, angular_spread_sim, theta)


#  diffa = numpy.zeros(60, complex)
#  diffpm = numpy.zeros(60, complex)
#  angle_sim = numpy.zeros(60, complex)
#  focal_ratio = numpy.zeros(60, complex)
#  
#  k = 0
#  for angular_spread_sim in numpy.arange(0.98,1.51,0.01):
#      (W, logW, logWfixed, n, tfdiff, tfexp, fscale, tscale, weight, Wnew, w, weight, avgdiffa, avgdiffpm)=process.extraction(N, M, dk, deltat, deltaf, reference_file, sample_file, thickness, convthickness, \
#                     setthetacutoff, processing, kstart, setkmin, paddingfactor, phasetolerance, setangularcentre, display, numRefl, polarisation_configuration, angular_spread_sim)
#      print k, angular_spread_sim,1/(2*numpy.tan(angular_spread_sim/2))
#      diffa[k]=avgdiffa
#      diffpm[k]=avgdiffpm
#      angle_sim[k]=angular_spread_sim
#      focal_ratio[k]=1/(2*numpy.tan(angular_spread_sim/2))
#      k+=1
# 
# def savedatafilespread(y, filename):
#     M=len(y)
#     fileobj=open(filename, 'w')
#     for i in range(0,M):
#         fileobj.write(str(angle_sim[i].real)+","+str(focal_ratio[i].real)+","+str(y[i])+"\n")
#     fileobj.close()     
#     
# 
# savedatafilespread(diffa.real, 'avgdar.txt')
# savedatafilespread(diffa.imag, 'avgdai.txt')
# savedatafilespread(diffpm.real, 'avgdpmr.txt')
# savedatafilespread(diffpm.imag, 'avgdpmi.txt')
        


P=0; Q=N/2

pylab.figure()
pylab.plot(fscale[P:Q], n[P:Q].real)
pylab.title('Refractive index')
pylab.grid(True)
pylab.legend()

pylab.figure()
pylab.plot(fscale[P:Q], n[P:Q].imag)
pylab.title('Extinction coefficient')
pylab.legend()
pylab.grid(True) 
 
pylab.figure()
pylab.title('Amplitude of transfer function')
pylab.plot(fscale[P:Q], logWfixed[P:Q].real)
pylab.plot(fscale[P:Q], weight[P:Q].real)

pylab.figure()
pylab.plot(fscale[P:Q], W[P:Q].real, label='old')
pylab.plot(fscale[P:Q], Wnew[P:Q].real, label='new')
pylab.title('Transfer function real')
pylab.legend()
pylab.grid(True) 

pylab.figure()
pylab.plot(fscale[P:Q], W[P:Q].imag, label='old')
pylab.plot(fscale[P:Q], Wnew[P:Q].imag, label='new')
pylab.title('Transfer function imag')
pylab.legend()
pylab.grid(True)

 
pylab.figure()
pylab.plot(fscale[P:Q], logWfixed[P:Q].imag)
pylab.plot(fscale[P:Q], logW[P:Q].imag)
pylab.title('Phase')
pylab.legend()
pylab.grid(True) 

pylab.figure()
pylab.plot(fscale[P:Q], w[P:Q].real)
pylab.title('Phase')
pylab.legend()
pylab.grid(True)

#pylab.figure()
#pylab.plot(tscale[0:len(x)], x[0:len(x)])
#pylab.plot(tscale[0:len(x)], y[0:len(x)])
#pylab.title('real')
#pylab.legend()
#pylab.grid(True)

pylab.figure()
pylab.plot(fscale[P:Q], tfdiff[P:Q].real)
pylab.plot(fscale[P:Q], tfdiff[P:Q].imag)
pylab.title('tfdiff')
pylab.legend()
pylab.grid(True)

pylab.figure()
pylab.plot(fscale[P:Q], tfexp[P:Q].real)
pylab.plot(fscale[P:Q], logWfixed[P:Q].real)
pylab.title('tfamplitude')
pylab.legend()
pylab.grid(True)

pylab.figure()
pylab.plot(fscale[P:Q], tfexp[P:Q].imag)
pylab.plot(fscale[P:Q], logWfixed[P:Q].imag)
pylab.title('tfphase')
pylab.legend()
pylab.grid(True)


#theta_max=2.0*models.theta_cutoff
#theta=0.5*theta_max*(1.0 + models.ix)
#profile=0.5*theta_max*models.angular_cutoff(theta, models.theta_cutoff, ops.freqk, models.kmin)


#theta=numpy.zeros(0.4/0.001)
#i=0
#for k in numpy.arange(0,0.4,0.001):
#    theta[i]=k
#    i+=1
    
#profile=models.angular_cutoffdonut(theta, models.theta_cutoff, ops.freqk, models.angularcentre)
#pylab.figure()
#pylab.plot(theta,profile)
#pylab.show()

#need to put no of r, into parameter



#ops.savedatafilefreq(n.real, 'ig03nR3sifp.txt', 1, deltaf)
#ops.savedatafilefreq(n.imag, 'eg03nR3sifp.txt', 1, deltaf)

1/0

namereffile=str(reference_file)
namereffile=namereffile.replace(".txt", "")
namereffile=namereffile.replace("Si/", "")

namesamfile=str(sample_file)
namesamfile=namesamfile.replace("Si/", "")
namesamfile=namesamfile.replace(".txt", "")

#ops.savedatafilefreq(n.real, 'i_nR'+str(int(1e6*thickness[0]))+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
#ops.savedatafilefreq(n.imag, 'e_nR'+str(int(1e6*thickness[0]))+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)


ops.savedatafilefreq(n.real, 'i_nRtpf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
ops.savedatafilefreq(n.imag, 'e_nRtpf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)

1/0

ops.savedatafilefreq(n.real, 'i2670noR_pf'+str(paddingfactor)+'.txt', 1, deltaf)
ops.savedatafilefreq(n.imag, 'e2670noR_pf'+str(paddingfactor)+'.txt', 1, deltaf)


1/0

namereffile=str(reference_file)
namereffile=namereffile.replace(".txt", "")
namereffile=namereffile.replace("Pb_data/", "")
namereffile=namereffile.replace("SiO2Pb00nm", "S")

namereffile=namereffile.replace("Reference", "R")
namesamfile=str(sample_file)
namesamfile=namesamfile.replace("Pb_data/", "")
namesamfile=namesamfile.replace(".txt", "")

#ops.savedatafilefreq(logW.imag, 'op_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
#ops.savedatafilefreq(logWfixed.imag, 'p_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
#ops.savedatafilefreq(logWfixed.real, 'a_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
ops.savedatafilefreq(n.real, 'i_pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
ops.savedatafilefreq(n.imag, 'e_pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)




1/0
namereffile=str(reference_file)
namereffile=namereffile.replace(".txt", "")
namereffile=namereffile.replace("Si_GaAs/", "")
namereffile=namereffile.replace("ref", "R")

namesamfile=str(sample_file)
namesamfile=namesamfile.replace("Si_GaAs/", "")
namesamfile=namesamfile.replace(".txt", "")
namesamfile=namesamfile.replace("gaase", "Em")
namesamfile=namesamfile.replace("gaasrec", "Re")
namesamfile=namesamfile.replace("gaaspara", "Pa")

#ops.savedatafilefreq(logW.imag, 'op_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
#ops.savedatafilefreq(logWfixed.imag, 'p_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
#ops.savedatafilefreq(logWfixed.real, 'a_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
ops.savedatafilefreq(n.real, 'i_pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
ops.savedatafilefreq(n.imag, 'e_pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)




1/0
namereffile=str(reference_file)
namereffile=namereffile.replace(".dat", "")
namereffile=namereffile.replace("HDPE_", "")
namereffile=namereffile.replace("Reference", "R")
namesamfile=str(sample_file)
namesamfile=namesamfile.replace("HDPE_after_filament", "AF")
namesamfile=namesamfile.replace("HDPE_focus", "F")
namesamfile=namesamfile.replace("HDPE_collimated", "C")
namesamfile=namesamfile.replace(".dat", "")
namesamfile=namesamfile.replace("HDPE_", "")

#ops.savedatafilefreq(logW.imag, 'op_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
#ops.savedatafilefreq(logWfixed.imag, 'p_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
#ops.savedatafilefreq(logWfixed.real, 'a_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
ops.savedatafilefreq(n.real, 'i_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)
ops.savedatafilefreq(n.imag, 'e_pt'+str(10*phasetolerance)+'pf'+str(paddingfactor)+'_'+namesamfile+'_'+namereffile+'.txt', 1, deltaf)


#1/0
#ops.savedatafilefreq(logW.imag, 'opnR'+namesamfile+'_'+namereffile+'.txt', paddingfactor, deltaf)
#ops.savedatafilefreq(logWfixed.imag, 'pnR'+namesamfile+'_'+namereffile+'.txt', paddingfactor, deltaf)
#ops.savedatafilefreq(logWfixed.real, 'anR'+namesamfile+'_'+namereffile+'.txt', paddingfactor, deltaf)
#ops.savedatafilefreq(n.real, 'inR'+namesamfile+'_'+namereffile+'.txt', paddingfactor, deltaf)
#ops.savedatafilefreq(n.imag, 'enR'+namesamfile+'_'+namereffile+'.txt', paddingfactor, deltaf)

#1/0

#file = open('fdbp'+namesamfile+'_'+namereffile+'.txt','w')
#for i in range(0, len(n)):
#    file.writelines(str(i*deltaf*1e-12)+","+str(fdbeamprofile[i].real)+"\n")
#file.close()

#1/0

#ops.savedatafiletime(y[i].real, 'y.txt', deltat)
#ops.savedatafiletime(x[i].real, 'x.txt', deltat)
    


