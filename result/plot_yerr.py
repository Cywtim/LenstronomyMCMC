# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 11:59:36 2022

@author: CHENG
"""

import numpy as np
import scipy as sp
from scipy import stats as st
from matplotlib import pyplot as plt

Nobs = 110
Nobs_LCDM = 1000
Nobs_wCDM = 150
Nobs_w0waCDM = 150
FlatLCDM_Nobs_list = np.array([10,30.50,300,500,1000])
LCDM_Nobs_list = np.array([10,30,70,90,110,200,300,400,500,600,700,900,1000])
wCDM_Nobs_list = np.array([10, 50, 100, 300, 1000])

is_LCDM = 0
is_wCDM = 2
is_w0waCDM = 2

component = 1
FlatLCDM_label_list = ["$\Omega_M$", "$w_0$", "$H_0$ km/s/Mpc"]
LCDM_label_list = ["$\Omega_M$", "$\Omega_\Lambda$", "$H_0$ km/s/Mpc"]
wCDM_label_list = ["$\Omega_M$", "$\Omega_\Lambda$", "$\omega_0$", "$H_0$ km/s/Mpc"]
FlatLCDM_truth_value = [0.3, -1, 72]
LCDM_truth_value = [0.3 , 0.7 , 72]
wCDM_truth_value = [0.3 , 0.7 , -1, 72]

mean_list = []
quantile_list = []
error_list = []
for N in range(len(FlatLCDM_Nobs_list)):
    try:
        chain = np.load("FlatwCDM_chain_%d.npy"%FlatLCDM_Nobs_list[N])
        #autocorr = np.load("wCDM_autocorr_%d.npy"%Nobs_list[N])
        confidence_level = 0.99
        mean = np.mean(chain,axis=0)
        print(mean)
        freedom = len(chain) - 1
        std = st.sem(chain,axis=0)
        confidence = st.t.interval(alpha=confidence_level, df=freedom,\
                               loc=mean, scale=std)
        quantile = np.quantile(chain,(0.16,0.84),axis=0)
        quantile_mean = (quantile - mean) 
        
    except:
        
        pass
    
    mean_list.append(mean)
    quantile_list.append(quantile)
    error_list.append(quantile - mean)

mean_list = np.array(mean_list)
quantile_list = np.array(np.array(quantile_list))
error_list = np.abs(error_list)
#fig = plt.figure(221)
#ax = plt.Subplot(fig)
plt.errorbar(wCDM_Nobs_list, mean_list[:,component],\
             yerr = [error_list[:,0,component],error_list[:,1,component]], \
                 linestyle=":",markersize=9,\
                     fmt="sb",capsize=3,\
                         ecolor="k",barsabove=True)
plt.fill_between(wCDM_Nobs_list, quantile_list[:,0,component],\
                 quantile_list[:,1,component],\
                     color="papayawhip",edgecolor='c',alpha=0.5)
plt.xlabel("Number of clusters", size=20)
plt.ylabel(FlatLCDM_label_list[component], size=20)
plt.semilogx(Nobs)
plt.ylim(55,85)
plt.axhline(FlatLCDM_truth_value[component],color='r',linestyle='-.',alpha=0.8)
plt.xticks(size=15)
plt.yticks(size=15)
