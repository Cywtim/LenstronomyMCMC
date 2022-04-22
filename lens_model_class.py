# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 18:29:53 2022

@author: CHENG
"""

import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
import astropy as ap
from astropy import cosmology,constants
import re


class lens_redshift_difference:
    
    def __init__(self,lens_model_list):
        '''
        This class is used to compare the difference bewteen different lens 
        models. It receives model names as input.
        
        Parameters
        ----------
        lens_model_list : array or list
        
            the name of lens model mentioned in 
            https://lenstronomy.readthedocs.io/en/latest/_modules/index.html

        Returns
        -------
        None.
        
        '''
        
        self.lens_model_list = lens_model_list
        
        self.lensModel = LensModel(lens_model_list=lens_model_list)
        
        
    def lens_solver(self,beta, kwargs_lens, search_window=80, min_distance=0.1):
        '''

        Parameters
        ----------
        beta : array or list
            The position of the source
        kwargs_lens : dict
            corresponding key words to the lens
            model given when define the class
        search_window : int, optional
            the search window of the solver. The default is 100.

        Returns
        -------
        theta_ra : float
            the RA of images
        theta_dec : float
            the DEC of images

        '''
        
        lensEquationSolver = LensEquationSolver(self.lensModel)
        
        theta_ra, theta_dec = lensEquationSolver.image_position_from_source(\
                    kwargs_lens=kwargs_lens, sourcePos_x=beta[0],sourcePos_y=beta[1],\
                        search_window=search_window, min_distance=0.1)
            
        return np.array([theta_ra, theta_dec])
    
    def lens_potential(self,beta, kwargs_lens, search_window=80, min_distance=0.1):
        '''

        Parameters
        ----------
        beta : array of list
            DESCRIPTION.
        kwargs_lens : dict
            corresponding key words to the lens
            model given when define the class
        search_window : int, optional
            the search window of the solver. The default is 100.

        Returns
        -------
        fermat_potential : array
            fermat potential for each image

        '''
        
        theta_x, theta_y = self.lens_solver(beta, kwargs_lens,search_window=search_window, min_distance=0.1)
            
        fermat_potential = self.lensModel.fermat_potential(x_image=theta_x,y_image=theta_y\
                    ,x_source=beta[0],y_source=beta[1],kwargs_lens=kwargs_lens)\
                    /((3600*180)**2)*np.pi**2
            
        return fermat_potential
    
    
    def lens_redshifts(self, beta, kwargs_lens, model_cosmology,\
                      zl, zs, pi_cosmology,search_window=80,min_distance=0.1):
        """
        

        Parameters
        ----------
        beta : list or ndarray
            The position of source
        kwargs_lens : dict
            lens parameters
        model_cosmology : string 
            model name (method name in astropy.cosmology)
        pi_cosmology : ndarray
            cosmological parameters in the form of
                        [zl, zs, other paramters]

        Returns
        -------
        redshifts : ndarray
            redshifts of each images

        """
                
        
        
        #l,zs,pi_cosmology = np.split(np.array(pi_cosmology),[1,2])
        
        #zl = pi_cosmology[0]
        #zs = pi_cosmology[1]
        #pi_cosmology = pi_cosmology[2:]
        
        pi_cosmology = str(pi_cosmology)[1:-1]
        pi_cosmology = pi_cosmology.split(',')
        pi_cosmology = ','.join(pi_cosmology)
        
        cosmo = eval('cosmology.' + model_cosmology+ '('+pi_cosmology+')' )
        Dl = cosmo.angular_diameter_distance(zl)
        Ds = cosmo.angular_diameter_distance(zs)
        Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
        
        fermat_potential_images = self.lens_potential(beta, kwargs_lens,\
                               search_window=search_window, min_distance=0.1)
        
        redshifts = (-(1+zl)*cosmo.H(zs)*Dl*Ds/Dls/constants.c).decompose().value\
                            * fermat_potential_images
        
        return np.array(redshifts)
    
    def lens_redshift_difference(self, beta, kwargs_lens, model_cosmology,\
                       zl,zs, pi_cosmology,search_window=80,min_distance=0.1):
        """
        

        Parameters
        ----------
        beta : list or ndarray
            The position of source
        kwargs_lens : dict
            lens parameters
        model_cosmology : string 
            model name (method name in astropy.cosmology)
        pi_cosmology : ndarray
            cosmological parameters in the form of
                        [zl, zs, other paramters]

        Returns
        -------
        redshift difference : ndarray
            redshift difference of each images

        """
                
        
        
        #l,zs,pi_cosmology = np.split(np.array(pi_cosmology),[1,2])
        
        #zl = pi_cosmology[0]
        #zs = pi_cosmology[1]
        #pi_cosmology = pi_cosmology[2:]
        
        pi_cosmology = str(pi_cosmology)[1:-1]
        pi_cosmology = pi_cosmology.split(',')
        pi_cosmology = ','.join(pi_cosmology)
        
        cosmo = eval('cosmology.' + model_cosmology+ '('+pi_cosmology+')' )
        Dl = cosmo.angular_diameter_distance(zl)
        Ds = cosmo.angular_diameter_distance(zs)
        Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
        
        fermat_potential_images = self.lens_potential(beta, kwargs_lens,\
                               search_window=search_window, min_distance=0.1)
        
        fermat_potential_images = fermat_potential_images.min() - \
                         fermat_potential_images.max()
        
        redshifts = (-(1+zl)*cosmo.H(zs)*Dl*Ds/Dls/constants.c).decompose().value\
                            * fermat_potential_images
        
        return np.array(redshifts)

        

        
        
        
        
        
        