import numpy as np
import matplotlib.pyplot as plt

#Quartz data 
data_quartz = np.loadtxt(r"C:\Users\celes\OneDrive\Documents\Documents\Uni Stuff\Fourth Year\Theo's Work\THz-TDS-main\Transfer Matrix Method\Matrix_methods\Crete measurements\Si_HRFZ_0.38mm_nitrogen_10min.txt")

#Chose this data set as they had the highest peak amplitude
t_s_raw = data_quartz[:,0]
A_s = data_quartz[:,1] 

#Thickness
d = 0.4e-3
