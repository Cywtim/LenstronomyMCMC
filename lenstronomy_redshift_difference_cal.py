# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 16:19:23 2022

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

#choice of lens mdoel
lens_nie = 'NIE'
lens_spep = 'SPEP'
lens_sie = 'SIE'
lens_sis = 'SIS'
lens_model_list = [lens_nie]
lens_model_set = [[lens_nie],[lens_sie], [lens_spep]]

#set up model
lensModel = LensModel(lens_model_list=lens_model_list)


#def parameter values of lens models
Nobs = 100
kwargs_nie = {'theta_E':40, 'e1':0.5, 'e2':-0.5,\
              's_scale':1, 'center_x':0.1, 'center_y':0.1}
kwargs_spep = {'theta_E': 40, 'e1': 0.5, 'e2': -0.5, 
               'gamma': 3, 'center_x': 0.1, 'center_y': 0.1}
kwargs_sie = {'theta_E': 40, 'e1': 0.5, 'e2': -0.5,\
                     'center_x': 0.1, 'center_y': 0.1}
kwargs_sis = {'theta_E':40,'center_x':0.1,'center_y':0.1}

kwargs_nie_stand = {'sigma_v':1500, 'e1':0.1, 'e2':-0.1,\
              's_scale':1, 'center_x':0.1, 'center_y':0.1}


sigma_v = np.random.randn(Nobs) * 150  + 1500
e1 = np.random.randn(Nobs)*0.1 + 0.1
e2 = -np.random.randn(Nobs)*0.1 - 0.1
s_scale = np.random.randn(Nobs)*0.1 + 1
center_x = np.random.randn(Nobs)*0.01 + 0.1
center_y = np.random.randn(Nobs)*0.01 + 0.1
kwargs_names = np.array(['sigma_v', 'e1', 'e2', 's_scale', 'center_x', 'center_y'])
lens_kwargs = np.array([sigma_v, e1, e2, s_scale, center_x, center_y]).T
zl_true = np.random.randn(Nobs)*0.1 + 1
zs_true = np.random.randn(Nobs)*0.4 + 4

Omega_M_true = 0.3
Omega_Lambda_true = 0.7
H0_true = 72

Delta_z_true = []
Delta_z_obs = []
beta = [0,0]
lens = lens_redshift_difference(lens_model_list)
cosmology_model = 'LambdaCDM'

kwargs_stand = [kwargs_nie_stand]
redshift_difference_stand = lens.NIE_redshifts(\
                beta, kwargs_stand, cosmology_model, 1, 4, \
                    [H0_true, Omega_M_true,Omega_Lambda_true],\
                  search_window=100,min_distance=3,solver='lenstronomy')
redshift_difference_stand = redshift_difference_stand.max() - \
                                     redshift_difference_stand.min()

for i in range(Nobs):
    
    lens_kwargs_list = [dict(zip(kwargs_names, lens_kwargs[i]))]
    
    temp = lens.NIE_redshifts(beta, lens_kwargs_list, cosmology_model,
                               zl_true[i], zs_true[i], [H0_true, Omega_M_true,
                                        Omega_Lambda_true],\
                                   search_window=100,min_distance=3,\
                                       solver="lenstronomy")
    #print(len(temp))
    temp = temp.max() - temp.min()
    
    Delta_z_true.append(temp)
    
    print(i,temp,",",)

Delta_z_true = np.array(Delta_z_true)

'''
start_time = time.time()

nie = lens_redshift_difference([lens_nie])

nie_r = nie.lens_redshifts([0.0,0.0],kwargs_lens_list,cosmo_name,cosmo_pi,\
                           search_window=80,min_distance=30)

end_time = time.time()

delta_time = end_time - start_time

print(nie_r,',',delta_time)
'''





