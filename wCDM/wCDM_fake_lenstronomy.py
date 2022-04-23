# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 23:43:24 2022

@author: CHENG
"""

import sys 
sys.path.append("..") 
import lens_model_class
from lens_model_class import lens_redshift_difference
import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver

eps_true = 0.2
Nobs = 10
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
kwargs_nie = {'theta_E':40, 'e1':0.5, 'e2':-0.5,\
              's_scale':0.1, 'center_x':0.1, 'center_y':0.1}
kwargs_spep = {'theta_E': .9, 'e1': 0.05, 'e2': 0.05, 
               'gamma': 2., 'center_x': 0.1, 'center_y': 0}
kwargs_shear = {'gamma1':0.5, 'gamma2':0.1, 'ra_0':0.0, 'dec_0':0.0}
kwargs_sis = {'theta_E':40,'center_x':0.5,'center_y':0.1}
kwargs_lens = [kwargs_nie]

zl_true = np.random.uniform(1, 4, size=Nobs)
zs_true = zl_true + np.random.uniform(2.5 , 3.5, size=Nobs)
Omega_M_true = 0.3
Omega_Lambda_true = 0.7
w0_true = -1
H0_true = 72

Delta_z_true = []
Delta_z_obs = []
lens = lens_redshift_difference(lens_model_list)
cosmology_model = 'wCDM'

for i in range(Nobs):
    temp = lens.lens_redshifts(beta, kwargs_lens, cosmology_model,\
            zl_true[i], zs_true[i], [H0_true, Omega_M_true, Omega_Lambda_true, w0_true])
    temp = temp.max() - temp.min()
    temp2 = temp  + err * np.random.normal(0, eps_true)
    
    Delta_z_true.append(temp)
    Delta_z_obs.append(temp2)
    
    print(i,",",temp,",",temp2)

#Delta_z_obs = Delta_z_true * ( 1 + ( np.random.normal(0,eps_true,size=Nobs) ) )
#np.random.random(Nobs)

np.save("..\data\wCDM10-8.npy",[zl_true,zs_true,Delta_z_true,Delta_z_obs])







