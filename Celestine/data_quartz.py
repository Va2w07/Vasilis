import numpy as np
import matplotlib.pyplot as plt

#Quartz data 
data_quartz = np.loadtxt(r"./data/Quartz_1mm_10min_nitrogen.txt")

#Chose this data set as they had the highest peak amplitude
t_s_raw = data_quartz[:,0]
A_s = data_quartz[:,1] 

#Thickness
d = 1e-3
