# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 12:55:52 2022

@author: CHENG
"""

import numpy as np
import scipy as sp
from scipy import stats as st

Nobs = 110
Nobs_LCDM = 100
Nobs_wCDM = 150
Nobs_w0waCDM = 150

is_LCDM = 0
is_wCDM = 2
is_w0waCDM = 2


try:
    LCDM_chain = np.load("LCDM_chain_%d.npy"%Nobs)
    LCDM_autocorr = np.load("LCDM_autocorr_%d.npy"%Nobs)
    
    wCDM_chain = np.load("wCDM_chain_%d.npy"%Nobs)
    wCDM_autocorr = np.load("wCDM_autocorr_%d.npy"%Nobs)
    
    w0waCDM_chain = np.load("w0waCDM_chain_%d.npy"%Nobs)
    w0waCDM_autocorr = np.load("w0waCDM_autocorr_%d.npy"%Nobs)

except:
    
    pass



if is_LCDM == 0:
    
    confidence_level_LCDM = 0.99
    mean_LCDM = np.mean(LCDM_chain,axis=0)
    freedom_LCDM = len(LCDM_chain) -1
    std_LCDM = st.sem(LCDM_chain,axis=0)
    confidence_LCDM = st.t.interval(alpha=confidence_level_LCDM, df=freedom_LCDM,\
                           loc=mean_LCDM, scale=std_LCDM)
    quantile_LCDM = np.quantile(LCDM_chain,(0.16,0.85),axis=0)
    
elif is_LCDM == 1:
    pass

#########################################################################
if is_wCDM == 0:
    
    confidence_level_wCDM = 0.99
    mean_wCDM = np.mean(wCDM_chain,axis=0)
    freedom_wCDM = len(wCDM_chain) -1
    std_wCDM = st.sem(wCDM_chain,axis=0)
    confidence_wCDM = st.t.interval(alpha=confidence_level_wCDM, df=freedom_wCDM,\
                           loc=mean_wCDM, scale=std_wCDM)
    quantile_wCDM = np.quantile(wCDM_chain,(0.16,0.85),axis=0)
    
elif is_wCDM == 1:
    pass

##########################################################################
if is_w0waCDM == 0:
    
    confidence_level_w0waCDM = 0.99
    mean_w0waCDM = np.mean(w0waCDM_chain,axis=0)
    freedom_w0waCDM = len(w0waCDM_chain) -1
    std_w0waCDM = st.sem(w0waCDM_chain,axis=0)
    confidence_w0waCDM = st.t.interval(alpha=confidence_level_w0waCDM,\
                        df=freedom_w0waCDM,loc=mean_w0waCDM, scale=std_w0waCDM)
    quantile_w0waCDM = np.quantile(w0waCDM_chain,(0.16,0.85),axis=0)
    
elif is_w0waCDM == 1:
    pass

