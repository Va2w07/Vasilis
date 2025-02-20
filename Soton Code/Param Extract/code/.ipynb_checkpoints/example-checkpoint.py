import sys
import ops
import numpy

def wFunc(maxi, length):
	c  = 100.0
	wVal = numpy.ones([length])
	for i in range(maxi,length):
		wVal[i] = numpy.exp(-(i - maxi)**2/(2*c**2))

	return wVal

#
###########################################################################
# Set data file names if appropriate
###########################################################################
reference_file = "data/2mmR.dat"
sample_file = "data/2mmS.dat"
#-------------------------------------------------------------------------
# Read the files, cols = number of columns in data file
#-------------------------------------------------------------------------
cols=2
# dt_ref will be set and null array returned for sd
input_scale=1.0
(dt_ref, x, ex)=ops.read_raw_datafile(cols, reference_file, input_scale)
(dt_data, y,ey)=ops.read_raw_datafile(cols, sample_file,    input_scale)
#-------------------------------------------------------------------------
# Reconfigure the raw data
#-------------------------------------------------------------------------
# set starting points of the data that are to be used in reference and 
# sample files
start_x=0
start_y=0
# set offsets in reference and data files if starting times not the same;
# signal is moved so that 'start' is 'offset'.
offset_x=0
offset_y=337
# set length of data to be processed, this should be power of 2 (no check!)
L=2**11
# chop the ends off the data, returns arrays of length L
x=ops.chop_shift_data(x,  start_x, offset_x, L)
y=ops.chop_shift_data(y,  start_y, offset_y, L)

xmaxi = numpy.argmax(x)
ymaxi = numpy.argmax(y)

x = x * wFunc(xmaxi, len(x))
y = y * wFunc(ymaxi, len(y))

#graph(x)
#graph(y)
#ey=ops.chop_shift_data(ey, start_y, offset_y, L)
#
#---------------------------------------------------------------------------
# Set sampling interval
#---------------------------------------------------------------------------
# set value derived from data file
#deltat=dt_ref
# if the file does not contain timing information deltat must be set here
# this is the value that will be used to generate synthetic data
# sampling interval 0.13ps, total duration 136.5ps
# Nyquist frequency 3.75THz, samples every ~3.7GHz
deltat=0.0194e-12

# interpolate to change sampling interval, deltat is the required value,
# normally only one of the following lines should be operative
#x=ops.interpolate(x, dt_ref,  deltat)
#y=ops.interpolate(y, dt_data, deltat)
#
# reduce the data length to reduce over-sampling
#newL=1024
#x= ops.decimate(x,  L, newL)
#y= ops.decimate(y,  L, newL)
#ey=ops.decimate(sd, L, newL)
#L=newL
#
#---------------------------------------------------------------------------
# Create synthetic data
#---------------------------------------------------------------------------
# synthetic reference
#L=1024
#x=models.simulate_reference(L, deltat)
#---------------------------------------------------------------------------
# Single layer parallel beam
#---------------------------------------------------------------------------
#layers=[(4.0-0.005J, 3.0e-3)]
#(T,y)=models.simulate_parallel(x, layers, deltat, 0.002)
#(T,y)=models.simulate_parallel(x, layers, deltat, 0.0)
# truncate to L points
#y=y[0:L]
# save files if wanted
#simulated_ref_file="simulated_ref.txt"
#simulated_data_file="simulated_data.txt"
#ops.savedatafile(2, deltat, x, x, simulated_ref_file,  1e-9)
#ops.savedatafile(2, deltat, y, y, simulated_data_file, 1e-9)
# set thickness used in processing
#thickness=array([layers[0][1]])
#---------------------------------------------------------------------------
# Sandwich, parallel beam
#---------------------------------------------------------------------------
#layers=[(2.0, 1.0e-3), (3.0-0.02J, 0.1e-3), (2.0, 1.0e-3)]
#(T,y)=models.simulate_parallel(x, layers, deltat, 0.002)
# truncate to L points
#y=y[0:L]
# save files
#simulated_ref_file="simulated_ref.txt"
#simulated_data_file="simulated_data.txt"
#ops.savedatafile(2, deltat, x, x, simulated_ref_file,  1e-9)
#ops.savedatafile(2, deltat, y, y, simulated_data_file, 1e-9)
#nbread=layers[0][0]
#thickness=array([layers[0][1], layers[1][1], layers[2][1]])
#---------------------------------------------------------------------------
# Converging beam
#---------------------------------------------------------------------------
#layers=[(2.0-0.005J, 2.0e-3)]
#theta_cutoff=0.2
#(T,y)=models.simulate_converging(x, 1.0, 2.0-0.005J, 2.0e-3, deltat, \
#                             theta_cutoff, 0.0, "gaussian")
#y=y[0:L]
#simulated_ref_file="simulated_ref.txt"
#simulated_data_file="simulated_data.txt"
#ops.savedatafile(2, deltat, x, x, simulated_ref_file,  1e-9)
#ops.savedatafile(2, deltat, y, y, simulated_data_file, 1e-9)
#thickness=array([layers[0][1]])
#-------------------------------------------------------------------------
# Set the errors to be used in the calculation.  These MUST be set
#-------------------------------------------------------------------------
# use values from file
#dy=array(ey)
# std dev from file may be too optimistic about errors
#dy=3.0*dy
# if std dev is not included in data file model the errors
# a simple model
print max(abs(y))
dy=(0.002*ones(L)*max(abs(y)) + 0.002*abs(y))
# a more sophisicated model
#dy=20.0*ops.noise_model_time(y)
# if nothing else set error constant
#dy=0.002*ones(L,float)
#-------------------------------------------------------------------------
# Set sample thickness used for processing, this is a list, even for one layer
#-------------------------------------------------------------------------
#thickness = [2.0e-3]
#thickness = [2.0e-2, 0.845e-2, 2.0e-2]

thickness = [2.0e-3]

#-------------------------------------------------------------------------
# Set sections of data to be ignored in the calculation
#-------------------------------------------------------------------------
# format is list of tuples consisting of start and finish data point numbers
ignore=[(0,0)]
#ignore=[(0,164), (1188,2047)]
#-------------------------------------------------------------------------
# Set various other parameters that may be needed in the chosen calculation
#-------------------------------------------------------------------------
# Set Lambda constant for impulse extending
#Lambda=1.0e-6
# Set refractive index of "bread" for GB sandwich
nquartz=2.1
#nbread=0.0
# Set cut-off angular scale for converging beam
#set_theta_cutoff(0.2,1e10)
#
#-------------------------------------------------------------------------
# Exit here if synthetic data is to be saved and not processed
#-------------------------------------------------------------------------
#sys.exit(0)
#-------------------------------------------------------------------------
# Set the processing operations for this calculation
#-------------------------------------------------------------------------
# Noise Model and Treatment:
#     constant noise no supression  1
#     conjugate gradient method     2
#
# Interpretation Model:
#     parallel single slab          1
#     parallel GB sandwich          2 (fits filling)
#     parallel filled roll          3 (fits outer layers)
#     converging beam               10
#
# Further processing
#     none                          0
#     display error bars on ri      1
#     fit ri to raw data            2
#     refine ri to smooth           3 (experimental)
#     refine ri to smooth           4 (new method, very experimental)

processing=(2, 1, 0)

# Choose what to see.
# raw data                                               1
# impulse response time domain                           2
# residuals in time                                      4
# residuals in frequency                                 8
# transfer function                                      16
# show unwrapped phase and log amplitude                 32
# show raw ri estimates                                  64
# debug; changes to ri                                   128
# show refined ri estimates                              256
# everything                                             511

display=511
display_errors=False
display_theory=False
