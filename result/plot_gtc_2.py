# -*- coding: utf-8 -*-
"""
Created on Tue May 24 23:22:00 2022

@author: CHENG
"""

import pygtc
import corner
import numpy as np
import matplotlib.pyplot as plt



def plot_gtc(Nobs, filename, universe_model="LCDM", plot_model=1):
    
    
    try:
        chain = np.load(filename)
    
    except:
        
        pass
    
    if universe_model == "LCDM":
        labels = ['$\Omega_M$', "$\Omega_\Lambda$", "$H_0$"]
        truths = [0.3, 0.7, 72]
    if universe_model == "FlatLCDM":
        labels = ['$\Omega_M$', "$H_0$"]
        truths = [0.3, 72]
    if universe_model == "wCDM":
        labels = ['$\Omega_M$', "$\Omega_\Lambda$", "w_0", "$H_0$"]
        truths = [0.3, 0.7, -1, 72]
    if universe_model == "FlatwCDM":
        labels = ['$\Omega_M$', "$w_0$", "$H_0$"]
        truths = [0.3, -1, 72]
    if universe_model == "w0waCDM":
        labels = ['$\Omega_M$', "$\Omega_\Lambda$", "$w_0$", "$w_a$", "$H_0$"]
        truths = [0.3, 0.7, -1, 0, 72]
    
    
    if plot_model == 0:
        corner.corner(chain[2000:],\
                    labels=labels,\
                        truths=truths,\
                            bins=40,hist_bin_factor=5,\
                                quantiles=(0.16, 0.84), levels=(1-np.exp(-0.5),),
                                max_n_ticks=20)
    elif plot_model == 1:
        pygtc.plotGTC(chains=[chain[2000:]],chainLabels=[universe_model]\
                      ,paramNames=labels,\
                          truths=truths,truthLabels='Standard',\
                              nContourLevels=3,figureSize='APJ_page',
                                  sigmaContourLevels=False,
                                  customTickFont={'size':12},
                                  customLabelFont={'size':15},
                                  customLegendFont={'size':20})
    else:
        pass

            
Nobs = 1000
filename = "FlatLCDM_chain_%d.npy"%Nobs
plot_gtc(Nobs, filename, universe_model="FlatLCDM")
