# -*- coding: utf-8 -*-
"""
Created on Tue May 24 23:22:00 2022

@author: CHENG
"""

import pygtc
import corner
import numpy as np
import matplotlib.pyplot as plt


Nobs_LCDM = 1000
Nobs_FlatLCDM = 10
Nobs_wCDM = 1000
Nobs_w0waCDM = 1000

try:
    FlatLCDM_chain = np.load("FlatLCDM_chain_%d.npy"%Nobs_FlatLCDM)
    
    FlatwCDM_chain = np.load("FlatwCDM_chain_%d.npy"%Nobs_FlatLCDM)
    
    LCDM_chain = np.load("LCDM_chain_%d-7.npy"%Nobs_LCDM)
    #LCDM_autocorr = np.load("LCDM_autocorr_%d.npy"%Nobs_LCDM)
    
    wCDM_chain = np.load("wCDM_chain_%d_2.npy"%Nobs_wCDM)
    #wCDM_autocorr = np.load("wCDM_autocorr_%d.npy"%Nobs_wCDM)
    
    w0waCDM_chain = np.load("w0waCDM_chain_%d_2.npy"%Nobs_w0waCDM)
    #w0waCDM_autocorr = np.load("w0waCDM_autocorr_%d.npy"%Nobs_w0waCDM)
    
    w0waCDM_sqrt_chain = np.load("w0waCDM_chain_sqrt_%d.npy"%Nobs_w0waCDM)
    w0waCDM__sqrt_autocorr = np.load("w0waCDM_autocorr_sqrt_%d.npy"%Nobs_w0waCDM)

except:
    
    pass


is_FlatLCDM = 2
is_FlatwCDM = 1
is_LCDM = 2
is_wCDM = 2
is_w0waCDM = 2


if is_FlatLCDM == 0:
    corner.corner(FlatLCDM_chain[2000:],\
                labels=['$\Omega_M$',"$H_0$"],\
                    truths=[0.3, 72],\
                        bins=40,hist_bin_factor=5,\
                            quantiles=(0.16, 0.84), levels=(1-np.exp(-0.5),),
                            max_n_ticks=20)
elif is_FlatLCDM == 1:
    pygtc.plotGTC(chains=[FlatLCDM_chain[2000:]],chainLabels=["wCDM"]\
                  ,paramNames=["$\Omega_M$", "$H_0$"],\
                      truths=(0.3, 72),truthLabels='Standard',\
                          nContourLevels=3,figureSize='APJ_page',
                              sigmaContourLevels=False,
                              customTickFont={'size':12},
                              customLabelFont={'size':15},
                              customLegendFont={'size':20})
else:
    pass

if is_LCDM == 0:
    corner.corner(LCDM_chain[2000:],\
                labels=['$\Omega_M$', '$\Omega_\Lambda$',"$H_0$"],\
                    truths=[0.3, 0.7, 72],\
                        bins=40,hist_bin_factor=5,\
                            quantiles=(0.16, 0.84), levels=(1-np.exp(-0.5),),
                            max_n_ticks=20)
elif is_LCDM == 1:
    pygtc.plotGTC(chains=[LCDM_chain[2000:]],chainLabels=["wCDM"]\
                  ,paramNames=["$\Omega_M$","$\Omega_\Lambda$", "$H_0$"],\
                      truths=(0.3,0.7, 72),truthLabels='Standard',\
                          nContourLevels=3,figureSize='APJ_page',
                              sigmaContourLevels=False,
                              customTickFont={'size':12},
                              customLabelFont={'size':15},
                              customLegendFont={'size':20})
else:
    pass

if is_FlatwCDM == 0:
    corner.corner(FlatwCDM_chain[2000:],\
                labels=['$\Omega_M$',"$\omega_0$", "$H_0$"],\
                    truths=[0.3, -1, 72],\
                        bins=40,hist_bin_factor=5,\
                            quantiles=(0.16, 0.84), levels=(1-np.exp(-0.5),),
                            max_n_ticks=20)
elif is_FlatwCDM == 1:
    pygtc.plotGTC(chains=[FlatwCDM_chain[2000:]],chainLabels=["wCDM"]\
                  ,paramNames=["$\Omega_M$", "$\omega_0$","$H_0$"],\
                      truths=(0.3, -1, 72),truthLabels='Standard',\
                          nContourLevels=3,figureSize='APJ_page',
                              sigmaContourLevels=False,
                              customTickFont={'size':12},
                              customLabelFont={'size':15},
                              customLegendFont={'size':20})
else:
    pass



if is_wCDM == 0:
    corner.corner(wCDM_chain[2000:],\
                labels=['$\Omega_M$', '$\Omega_\Lambda$',"$\omega$","$H_0$"],\
                    truths=[0.3, 0.7, -1, 72],\
                        bins=40,hist_bin_factor=5,\
                            quantiles=(0.16, 0.84), levels=(1-np.exp(-0.5),))
elif is_wCDM == 1:
    pygtc.plotGTC(chains=[wCDM_chain[2000:]],chainLabels=["wCDM"]\
                  ,paramNames=["$\Omega_M$","$\Omega_\Lambda$",\
                               "$\omega_0$", "$H_0$"],\
                      truths=(0.3,0.7,-1, 72),truthLabels='Standard',\
                          nContourLevels=3,figureSize='APJ_page',
                          customTickFont={'size':12},
                          customLabelFont={'size':15},
                          customLegendFont={'size':20})
else:
    pass

if is_w0waCDM == 0:
    corner.corner(w0waCDM_chain[2000:],\
                labels=['$\Omega_M$', '$\Omega_\Lambda$', "$w_0$", "$w_a$", "$H_0$"],\
                    truths=[0.3, 0.7, -1, 0, 72],\
                        bins=40,hist_bin_factor=5,\
                            quantiles=(0.16, 0.84), levels=(1-np.exp(-0.5),))
elif is_w0waCDM == 1:
    pygtc.plotGTC(chains=[w0waCDM_chain[2000:]],chainLabels=["w0waCDM"]\
                  ,paramNames=['$\Omega_M$', '$\Omega_\Lambda$', "$w_0$", "$w_a$", "$H_0$"],\
                      truths=(0.3, 0.7, -1, 0, 72),truthLabels='Standard',\
                          nContourLevels=3,figureSize='APJ_page',
                          customTickFont={'size':12},
                          customLabelFont={'size':15},
                          customLegendFont={'size':20})
else:
    pass