import numpy as np
import matplotlib.pyplot as plt

#Reference data
data_ref = np.loadtxt(r"./data/Reference_15min_nitrogen (1).txt")    

#Quartzsi data 
data_quartzsi = np.loadtxt(r"./data/Quartz_Si_15min_nitrogen.txt")


t_r_raw = data_ref[:,0]
A_r = data_ref[:,1]

t_s_raw = data_quartzsi[:,0]
A_s = data_quartzsi[:,1] 

#Thickness
d1 = 1e-3
d2 = 0.4e-3
