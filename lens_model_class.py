# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 18:29:53 2022

@author: CHENG

"""

import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from astropy import cosmology,constants
import math


class lens_redshift_difference:
    
    def __init__(self,lens_model_list):
        '''
        This class is used to compare the difference between different lens
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
        
        self.lens_model_list = lens_model_list #The lens model's name
        
        self.lensModel = LensModel(lens_model_list=lens_model_list) #build up the defined lens
        
        
    def lens_solver(self,beta, kwargs_lens, search_window=80,
                    min_distance=0.1,solver='lenstronomy'):
        '''

        This function is to solve the image position

        Parameters
        ----------
        beta : array or list
            The position of the source
        kwargs_lens : dict
            corresponding keywords to the lens
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
                        search_window=search_window,min_distance=min_distance,\
                            solver=solver)
            
        #print([theta_ra, theta_dec])
        return np.array([theta_ra, theta_dec])
    
    def lens_potential(self,beta, kwargs_lens, search_window=80,
                       min_distance=0.1,solver='lenstronomy'):
        '''

        Compute the potential for each image

        Parameters
        ----------
        beta : array of list
            DESCRIPTION.
        kwargs_lens : dict
            corresponding keywords to the lens
            model given when define the class
        search_window : int, optional
            the search window of the solver. The default is 100.

        Returns
        -------
        fermat_potential : array
            fermat potential for each image

        '''
        
        theta_x, theta_y = self.lens_solver(beta, kwargs_lens,
                                            search_window=search_window, 
                                            min_distance=min_distance,\
                                                solver=solver)
            
        fermat_potential = self.lensModel.fermat_potential(x_image=theta_x,
                     y_image=theta_y,x_source=beta[0],y_source=beta[1],
                     kwargs_lens=kwargs_lens)/((3600*180)**2)*np.pi**2
            
        return fermat_potential
    
    
    def lens_redshifts(self, beta, kwargs_lens, model_cosmology,\
                      zl, zs, pi_cosmology,
                      search_window=80,min_distance=0.1,solver='lenstronomy'):
        """

        Calculate the redshift in each image

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
                        [zl, zs, other parameters]

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
                               search_window=search_window,\
                                   min_distance=min_distance,solver=solver)
        
        redshifts = (-(1+zl)*cosmo.H(zs)*Dl*Ds/Dls/constants.c).decompose().value\
                            * fermat_potential_images
        
        return np.array(redshifts)
    
    def NIE_timedelay(self, beta, kwargs_lens, model_cosmology,\
                      zl, zs, pi_cosmology,
                      search_window=80,min_distance=0.1,solver='lenstronomy'):
        """

        The time delay in each image for non-singularity Isotropy Ellipse (NIE) lens model

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
                        [zl, zs, other parameters]

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
        
        c = constants.c.value / 1000
        kwargs_lens[0]['theta_E'] = math.degrees(4*np.pi*(kwargs_lens[0]['sigma_v']/c)**2 * (Dls/Ds).value) * 3600
        kwargs_lens[0].pop('sigma_v')
        fermat_potential_images = self.lens_potential(beta, kwargs_lens,\
                               search_window=search_window,\
                                   min_distance=min_distance,solver=solver)
        
        redshifts = (-(1+zl)*Dl*Ds/Dls/constants.c).decompose().value\
                            * fermat_potential_images / 60 / 60 / 24
        
        return np.array(redshifts)    
    

    
    def NIE_redshifts(self, beta, kwargs_lens, model_cosmology,\
                      zl, zs, pi_cosmology,
                      search_window=80,min_distance=0.1,solver='lenstronomy'):
        """

        The redshift in each image for non-singularity Isotropy Ellipse (NIE) lens model

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
                        [zl, zs, other parameters]

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
        
        c = constants.c.value / 1000
        kwargs_lens[0]['theta_E'] = math.degrees(4*np.pi*(kwargs_lens[0]['sigma_v']/c)**2 * (Dls/Ds).value) * 3600
        kwargs_lens[0].pop('sigma_v')
        fermat_potential_images = self.lens_potential(beta, kwargs_lens,\
                               search_window=search_window,\
                                   min_distance=min_distance,solver=solver)
        
        redshifts = (-(1+zl)*cosmo.H(zs)*Dl*Ds/Dls/constants.c).decompose().value\
                            * fermat_potential_images
        
        return np.array(redshifts)    
    
    def lens_redshift_difference(self, beta, kwargs_lens, model_cosmology,\
                       zl,zs, pi_cosmology,\
                       search_window=80,min_distance=0.1,solver='lenstronomy'):
        """
        !May not work!
        The biggest redshift difference between images

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
                        [zl, zs, other parameters]

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
                               search_window=search_window,\
                                   min_distance=min_distance,\
                                       solver=solver)
        
        fermat_potential_images = fermat_potential_images.min() - \
                         fermat_potential_images.max()
        
        redshifts = (-(1+zl)*cosmo.H(zs)*Dl*Ds/Dls/constants.c).decompose().value\
                            * fermat_potential_images
        
        return np.array(redshifts)

        

        
        
        
        
        
        