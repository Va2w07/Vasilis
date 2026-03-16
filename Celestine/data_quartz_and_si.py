import numpy as np
import matplotlib.pyplot as plt

#Quartz data 
data_quartz1 = np.loadtxt(r"./data/QREF_12.9um_nitrogen_3min.txt")
data_quartz2 = np.loadtxt(r"./data/QREF_12.9um_nitrogen_7min.txt")
data_quartz3 = np.loadtxt(r"./data/QREF_12.9um_nitrogen_10min.txt")    
data_quartz4 = np.loadtxt(r"./data/QREF_12.9um_start.txt")   

#Chose this data set as they had the highest peak amplitude
t_s_raw = data_quartz3[:,0]
A_s = data_quartz3[:,1] 

#Thickness
d_q = 1e-3+119.2e-6+1e-3 #Thicknesses of each induvidual layer

#Reference data
data_ref1 = np.loadtxt(r"./data/reference_3min nitrogen.txt")
data_ref2 = np.loadtxt(r"./data/reference_10min_nitrogen_2.txt")
data_ref3 = np.loadtxt(r"./data/reference_10min_nitrogen_3.txt")
data_ref4 = np.loadtxt(r"./data/reference_10min_nitrogen.txt")
data_ref5 = np.loadtxt(r"./data/reference_start.txt")


t_r_raw = data_ref3[:,0]
A_r = data_ref3[:,1]


#Silicon data
data_si1 = np.loadtxt(r"./data/Si_HRFZ_0.38mm_nitrogen_3min.txt")
data_si2 = np.loadtxt(r"./data/Si_HRFZ_0.38mm_nitrogen_6min.txt")
data_si3 = np.loadtxt(r"./data/Si_HRFZ_0.38mm_nitrogen_10min.txt")
data_si4 = np.loadtxt(r"./data/Si_HRFZ_0.38mm_nitrogen_start.txt")


t_si_raw = data_si3[:,0]
A_si = data_si3[:,1]

#Thickness
d_si = 380e-6
