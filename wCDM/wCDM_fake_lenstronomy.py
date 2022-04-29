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
lens = lens_redshift_difference(lens_model_list)

#generate several lens parameters
kwargs = np.array(["theta_e","e1","e2","s_scale","center_x","center_y"])
theta_E = np.random.randn(Nobs) * 5  + 40
e1 = np.random.randn(Nobs)*0.1 + 0.5
e2 = -np.random.randn(Nobs)*0.1 - 0.5
s_scale = np.random.randn(Nobs)*0.1 + 1
center_x = np.random.randn(Nobs)
center_y = np.random.randn(Nobs)
kwargs_names = np.array(['theta_E', 'e1', 'e2', 's_scale', 'center_x', 'center_y'])
lens_kwargs = np.array([theta_E, e1, e2, s_scale, center_x, center_y]).T



#form of the lens configuration
'''
kwargs_nie = {'theta_E':40, 'e1':0.5, 'e2':-0.5,\
              's_scale':0.1, 'center_x':0.1, 'center_y':0.1}
kwargs_lens = [kwargs_nie]
'''

zl_true = np.random.uniform(1, 4, size=Nobs)
zs_true = zl_true + np.random.uniform(2.5 , 3.5, size=Nobs)
Omega_M_true = 0.3
Omega_Lambda_true = 0.7
w0_true = -1
H0_true = 72

Delta_z_true = []
Delta_z_obs = []

cosmology_model = 'wCDM'

for i in range(Nobs):
    
    lens_kwargs_list = [dict(zip(kwargs_names, lens_kwargs[i]))]
    
    temp = lens.lens_redshifts(beta, lens_kwargs_list, cosmology_model,\
            zl_true[i], zs_true[i], [H0_true, Omega_M_true, Omega_Lambda_true, w0_true])
    temp = temp.max() - temp.min()
    temp2 = temp  + err * np.random.normal(0, eps_true)
    
    Delta_z_true.append(temp)
    Delta_z_obs.append(temp2)
    
    print(i,",",temp,",",temp2)

#Delta_z_obs = Delta_z_true * ( 1 + ( np.random.normal(0,eps_true,size=Nobs) ) )
#np.random.random(Nobs)

np.save("..\data\wCDM%d-8.npy"%Nobs,[zl_true,zs_true,Delta_z_true,Delta_z_obs])
np.save("..\data\wCDM%d-8_lens_parameters.npy"%Nobs,lens_kwargs)







