# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 23:43:24 2022

@author: CHENG
"""
import sys 
sys.path.append("..") 
import lens_model_class
import lens_model_class
from lens_model_class import lens_redshift_difference
import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver

eps_true = 1
Nobs = 1000
beta = [0,0]
err = 1e-8

#choice of lens mdoel
lens_nie = 'NIE'
lens_spep = 'SPEP'
lens_shear = 'SHEAR'
lens_sis = 'SIS'
lens_model_list = [lens_nie]

#set up model
lensModel = LensModel(lens_model_list=lens_model_list)


#def parameter values of lens models

#generate several lens parameters
kwargs = np.array(["theta_e","e1","e2","s_scale","center_x","center_y"])
theta_E = np.random.randn(Nobs) * 10  + 100
sigma_v = np.random.randn(Nobs) * 150  + 1500
e1 = np.random.randn(Nobs)*0.01 + 0.1
e2 = -np.random.randn(Nobs)*0.01 - 0.1
s_scale = np.random.randn(Nobs)*0.1 + 1
center_x = np.random.randn(Nobs)*0.1
center_y = np.random.randn(Nobs)*0.1
kwargs_names = np.array(['sigma_v', 'e1', 'e2', 's_scale', 'center_x', 'center_y'])
lens_kwargs = np.array([sigma_v, e1, e2, s_scale, center_x, center_y]).T

zl_true = np.random.uniform(1, 2, size=Nobs)
zs_true = zl_true + np.random.uniform(2 , 4, size=Nobs)
Omega_M_true = 0.3
Omega_Lambda_true = 1 - Omega_M_true
w0 = -1
H0_true = 72

Delta_z_true = []
Delta_z_obs = []
lens = lens_redshift_difference(lens_model_list)
cosmology_model = 'wCDM'
ratio_list = []

for i in range(Nobs):
    
    lens_kwargs_list = [dict(zip(kwargs_names, lens_kwargs[i]))]
    
    temp = lens.NIE_redshifts(beta, lens_kwargs_list, cosmology_model,
                               zl_true[i], zs_true[i], [H0_true, Omega_M_true,
                                        1-Omega_M_true, w0],\
                                   search_window=100,min_distance=3,\
                                       solver="lenstronomy")
    #print(len(temp))
    temp = temp.max() - temp.min(-1)
    temp2 = temp  + err * np.random.normal(0, eps_true)
    
    Delta_z_true.append(temp)
    Delta_z_obs.append(temp2)

    ratio_list.append(temp/temp2)
    print(i,temp,",",temp2,",",temp/temp2)

np.save("..\data\FlatwCDM%d-8.npy"%Nobs,[zl_true,zs_true,Delta_z_true,Delta_z_obs])
np.save("..\data\FlatwCDM%d-8_lens_parameters.npy"%Nobs,lens_kwargs)





