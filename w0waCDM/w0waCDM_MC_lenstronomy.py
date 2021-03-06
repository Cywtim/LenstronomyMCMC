# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:16:05 2022

@author: CHENG
"""
import sys 
sys.path.append("..") 
import lens_model_class
import lens_model_class
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as opt
import h5py
import emcee
import corner
import lens_model_class
from lens_model_class import lens_redshift_difference
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver

global number
number = 0

def lnprior(p):
    # The parameters are stored as a vector of values, so unpack them
    # If the parameters are out of expected range, set the value -inf
    # If not return 0 and calculate likehood
    Omega_M, Omega_Lambda, w0, wa, H0  = p
    # We're using only uniform priors, and only eps has a lower bound
    if Omega_M < 0.0  or Omega_M > 1.0 :
        return - np.inf
    elif Omega_Lambda < 0.0 or Omega_Lambda > 1.0:
        return - np.inf
    elif w0 < -2.0 or w0 > 0.0:
        return - np.inf
    elif wa < -1.0 or wa > 1.0:
        return - np.inf
    elif H0 < 50 or H0 > 90:
        return - np.inf

    return 0

def lnlike(p, zl, zs, Deltaz, z_err):
    #likehood function, use mean square function here
    
    
    global number
    
    
    Omega_M, Omega_Lambda, w0, wa, H0 = p
    
    model = np.array([])
    for i in range(len(zl)):
        
        lens_kwargs_list = [dict(zip(kwargs_names, lens_kwargs[i]))]
        
        temp = lens.NIE_redshifts(beta, lens_kwargs_list, cosmology_model,\
                zl_true[i], zs_true[i], [H0, Omega_M, Omega_Lambda, w0, wa],\
                    search_window=100,min_distance=3,solver="lenstronomy")
        temp = temp.max() - temp.min()
        
        model = np.append(model,temp)
        
    
    # the likelihood is sum of the lot of normal distributions
    
    denom = np.power(z_err,2)
    lp = - 0.5 * sum( np.power((Deltaz - model),2)/denom + np.log(denom))
    
    #number = number + 1
    #if number%100 == 0:
    #    print("i:likehood",number,",",p,",lp=",lp)
    
    return lp



def lnprob(p, zl_true, zs_true, Delta_z_obs, z_err):
    
    global number
    
    lp = lnprior(p) # test whether the vector is not out range
    
    if not np.isfinite(lp):
        
        #number = number + 1 
        #if number% 100 == 0:
        #    print("i:prior",number,",",p,",lp=-inf")    
        return - np.inf
    
    return lp + lnlike(p, zl_true, zs_true, Delta_z_obs, z_err) # If not compute likehood




zl_true,zs_true,Delta_z_true,Delta_z_obs = np.load("../data/w0waCDM10-8.npy")
lens_kwargs = np.load("../data/w0waCDM10-8_lens_parameters.npy")
kwargs_names = np.array(['sigma_v', 'e1', 'e2', 's_scale', 'center_x', 'center_y'])
z_err = 1e-8
Omega_M_true = 0.3
Omega_Lambda_true = 0.7
H0_true = 72
w0_true = -1.0
wa_true = 0.0
beta = [0, 0]

#nll = lambda *args: -lnlike(*args)
#result = opt.minimize(nll, [Omega_M_true, Omega_Lambda_true, H0_true, eps_true],
#                      args=(zl_true,zs_true,Delta_z_obs))

nwalkers, ndim = 10, 5

# create a nearby proir 
p0 = np.array([[Omega_M_true, Omega_Lambda_true, w0_true, wa_true, H0_true ]])
p1 = np.array([[0.0, 0.0, 0.0, 0.01, 0.0 ]])
p0 = np.array([ (p0+p1)*( 1 + 1.e-9*np.random.randn(ndim))\
               for i in range(nwalkers)]).reshape(nwalkers,ndim)


filename = 'wCDM10.h5'
backend = emcee.backends.HDFBackend(filename)
backend.reset(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,\
                                args=[zl_true, zs_true, Delta_z_obs, z_err],\
                                    backend=backend)
    
#track how the average autocorrelation time estimate changes
steps = 10000
index = 0
autocorr = np.empty(steps)

# This will be useful to testing convergence
old_tau = np.inf

#choice of lens mdoel
lens_nie = 'NIE'
lens_model_list = [lens_nie]

#def parameter values of lens models
'''
kwargs_nie = {'theta_E':40, 'e1':0.5, 'e2':-0.5,\
              's_scale':0.1, 'center_x':0.1, 'center_y':0.1}
kwargs_lens = [kwargs_nie]
'''
cosmology_model = 'wCDM'

lens = lens_redshift_difference(lens_model_list)


#sample for up to max_n steps
for sample in sampler.sample(p0, iterations=steps, progress=True, tune=True):
    # Only check convergence every 100 steps
    if sampler.iteration % 100:
        continue
    
    # Compute the autocorrelation time so far
    # Using tol=0 means
    #always get an estimate even if it isn't trustworthy
    tau = sampler.get_autocorr_time(tol=0)
    autocorr[index] = np.mean(tau)
    #print(tau)
    index = index + 1
    
    # Check convergence
    converged = np.all(tau * 100 < sampler.iteration)
    converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
    if converged:
        np.save("../result/w0waCDM_chain.npy", sampler.flatchain)
        np.save("../result/w0waCDM_autocorr.npy",autocorr)
        break
    old_tau = tau
        

np.save("../result/w0waCDM_chain.npy", sampler.flatchain)
np.save("../result/w0waCDM_autocorr.npy", autocorr)  
    

'''
pos,prob,state = sampler.run_mcmc(p0 , 5000, progress=True,\
                                  tune=True, store=True)

result = sampler.flatchain[np.where(np.sum(sampler.flatchain,axis=1)<1)\
                           or np.where(sampler.flatchain.any>1)\
                               or np.where(sampler.flatchain.any<0)]
pygtc.plotGTC(chains=[result[2000:]],chainLabels=["wCDM"]\
              ,paramNames=["$\Omega_M$","$\Omega_\Lambda$","$\omega_0$"],\
                  truths=(0.3,0.7,-1),truthLabels='Standard',\
                      nContourLevels=3,figureSize='APJ_page')
    
tmp = corner.corner(result[2000:],\
            labels=['$\Omega_M$', '$\Omega_\Lambda$',"$\omega$"],\
                truths=[Omega_M_true, Omega_Lambda_true, -1],\
                    bins=40,hist_bin_factor=5,\
                        quantiles=(0.16, 0.84))


'''


'''
tmp = pygtc.plotGTC(result,sigmaContourLevels=True,\
                    truths=[0.3,0.7],\
                        paramNames=['$\Omega_M$','$\Omega_\Lambda$'],\
                            figureSize='APJ_page',nBins=50,nContourLevels=3)
'''



'''
plt.figure(2)
plt.plot(autocorr[autocorr>0])
plt.xlabel('steps')
plt.ylabel('autocorrelation time $\tau$')
'''


