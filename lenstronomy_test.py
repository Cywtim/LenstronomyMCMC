# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 16:45:39 2022

@author: CHENG
"""

from lens_model_class import lens_redshift_difference
import numpy as np
import matplotlib.pyplot as plt
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
import time
import copy

#choice of lens mdoel
lens_nie = 'NIE'

lens_model_list = [lens_nie]

#set up model
lensModel = LensModel(lens_model_list=lens_model_list)


#def parameter values of lens models
kwargs_nie1 = {'theta_E':40, 'e1':0.11, 'e2':-0.09,\
              's_scale':1, 'center_x':0, 'center_y':0}
kwargs_nie2 = {'sigma_v':800, 'e1':0.1, 'e2':-0.2,\
              's_scale':0.1, 'center_x':0.10, 'center_y':0.5}



kwargs_lens_list = np.array([kwargs_nie2])



beta = [0.0, 0.0]

cosmo_name = 'w0waCDM'
zl = np.array([0.68])
zs = np.array([1.74])
cosmo_pi = [72,0.3,0.7,-1,0]
s_time = time.time()
lens_model = lens_redshift_difference(lens_model_list)

"""
in_time = time.time()

timedelay_list = lens_model.NIE_timedelay(\
                beta, kwargs_lens_list, cosmo_name, zl, zs, cosmo_pi,\
                  search_window=100,min_distance=3,solver='lenstronomy')
e_time = time.time()

print(timedelay_list.max() - timedelay_list.min())
"""


#print(redshift_difference,',',e_time - s_time)
#redshift_difference = redshift_difference.max() - redshift_difference.min()

wa_list = np.linspace(-1, 1)
redshift_difference_list = []
for i in range(len(wa_list)):
    cosmo_pi[-1] = wa_list[i]
    kwargs_lens = copy.deepcopy(kwargs_lens_list)
    redshift_difference = lens_model.NIE_redshifts(\
                    beta, kwargs_lens, cosmo_name, 0.5, 4, cosmo_pi,\
                      search_window=100,min_distance=3,solver='lenstronomy')
    redshift_difference = redshift_difference.max() - redshift_difference.min()
    redshift_difference_list.append(redshift_difference)


cosmo_pi[-1] = 0
kwargs_lens = copy.deepcopy(kwargs_lens_list)
redshift_difference = lens_model.NIE_redshifts(\
                beta, kwargs_lens, cosmo_name, 0.5, 4, cosmo_pi,\
                  search_window=100,min_distance=3,solver='lenstronomy')
standard = redshift_difference.max() - redshift_difference.min()

plt.figure()
plt.plot(wa_list, np.abs(redshift_difference_list-standard))
plt.xlabel("$w_a$")
plt.ylabel("Redshift difference")



"""
rsd = []
arg_range = np.array([60,130])
for arg in arg_range:
    kwargs_nie1['theta_E'] = arg
    redshift_difference_list = lens_model.lens_redshifts(\
                      beta, kwargs_lens_list, cosmo_name, 1, 4, cosmo_pi,\
                          search_window=300,solver='lenstronomy')
    
    rsd.append(redshift_difference_list.max() - redshift_difference_list.min())
"""

