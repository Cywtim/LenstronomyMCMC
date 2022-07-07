# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 23:55:31 2022

@author: CHENG
"""

import numpy as np
import scipy as sp
import astropy as ap
from scipy.optimize import root
from astropy import cosmology,constants
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_bvp, solve_ivp

class redshift_difference_func:
    """
    
    """
    
    def __init__(self,sigma=4.024922359499621,theta_lens=2.62,\
                 s=1,ellipticity=0.14):
        """
        Parameters
        ----------
        beta : 2d array, optional
            image position. 
        theta_lens : float or array, optional
            lens ratation. 
        sigma : float, optional.
            velocity dispersion of lens.
        s : float, optional
            core size of lens. 
        ellipticity : float, optional
            ellipticity. 
        
        """
        self.sigma = sigma
        
        self.theta_lens = theta_lens
        
        self.s = s
        
        self.ellipticity = ellipticity        
        
        
    def position_trans(self,theta,beta):
        
        
        rotation = [[np.cos(self.theta_lens), -np.sin(self.theta_lens)],
                    [np.sin(self.theta_lens), np.cos(self.theta_lens)]]
        
        if len(theta.shape) <= 2: 
            return np.matmul(rotation, np.array([ theta[0]-beta[0], theta[1]-beta[1]]))
        elif len(theta.shape) == 3:
            return np.array([ np.matmul(rotation, theta[:,0][i] - beta)\
                    for i in range(len(theta[:,0]))])
    
    
    def potential_tau(self,theta,beta,to_arc=0):
        
        x_,y_ = self.position_trans(theta,beta)
        tau = self.b*np.sqrt(self.s**2+(1-self.ellipticity)*x_**2+
                              (1+self.ellipticity)*y_**2)
        
        if to_arc == 0:
            return tau
        else:
            return tau/(3600*180/np.pi)**2 
    
    
    def d_potential_tau(self,beta,theta):
        
        
        x_,y_ = self.position_trans(theta, beta)
        d_tau = np.array([ beta[0]-theta[0]+self.b**2/self.potential_tau(theta)*\
            ((1-self.ellipticity)*np.cos(self.theta_lens)*x_+\
             (1+self.ellipticity)*np.sin(self.theta_lens)*y_)  ,\
            beta[1]-theta[1]+self.b**2/self.potential_tau(theta)*\
            ((1-self.ellipticity)*(-np.sin(self.theta_lens))*x_+\
             (1+self.ellipticity)*np.cos(self.theta_lens)*y_)  ])
            
        return d_tau
        
        
    def solve_images(self,Dls,Ds,beta,err=1e-2,radi=5,grid=3,number_output=False):
        """
        This function calculate the image positions 
        for a given source position beta
        
        input: source position
        output: image position list
        
        """        
        
        b = (4*Dls/Ds*self.sigma**2).decompose().value
        if type(b) != np.ndarray:
            b = np.array([b])
        #print(self.b)
        
        def equation(theta):   
            '''
            the gradiant of time delay surface/
            lens eqation/
            Eqn. (6)
            '''

            
            x_,y_ = self.position_trans(theta,beta)
            return np.array([ beta[0]-theta[0]+self.b**2/self.potential_tau(theta,beta)*\
                    ((1-self.ellipticity)*np.cos(self.theta_lens)*x_+\
                     (1+self.ellipticity)*np.sin(self.theta_lens)*y_)  ,\
                    beta[1]-theta[1]+self.b**2/self.potential_tau(theta,beta)*\
                    ((1-self.ellipticity)*(-np.sin(self.theta_lens))*x_+\
                     (1+self.ellipticity)*np.cos(self.theta_lens)*y_)  ])
        
        
        
        image_list = []
        for k in range(len(b)):
            
            theta_list = np.array([[]])
            self.b = b[k]
            
            for i in np.linspace(beta[0]-radi,beta[1]+radi,grid):
                for j in np.linspace(beta[0]-radi,beta[1]+radi,grid): 
                    '''
                    search in the lens plane
                    '''
                
                    
                    temp_theta = root(equation,np.array([i,j])).x
                    '''
                    solve the closest root near [i,j]
                    '''
    
                    '''
                    if the root is new, append the root in the list
                    '''                
                    if np.abs(equation(temp_theta)[0]) <= 1e-10 and np.abs(equation(temp_theta)[1]) <= 1e-10:
                        if len(theta_list[0]) == 0:
                            theta_list = np.append(theta_list,temp_theta)
                            theta_list = theta_list.reshape(1,2)
                        elif ((np.abs(temp_theta[0] - theta_list[:,0])>err) | \
                            (np.abs(temp_theta[1] - theta_list[:,1])>err)).all():
                            theta_list = np.append(theta_list,[temp_theta],axis=0)
                            
                

            '''
            sort the image list
            '''
            for i in range(len(theta_list)):
                
                for j in range(len(theta_list)-i-1):
                    
                    if theta_list[j][0] > theta_list[j+1][0]:
                        theta_list[[j,j+1],:] = theta_list[[j+1,j],:];
                    elif theta_list[j][0] == theta_list[j+1][0] and theta_list[j][1] > theta_list[j+1][1]:
                        theta_list[[j,j+1],:] = theta_list[[j+1,j],:];
            
            image_list.append(theta_list)
            
        if len(image_list) == 1:
            return np.array(image_list[0])
        else:
            return np.array(image_list)
        
        
    
    
    def LCDM(self,zl,zs,beta,ith=[0,1,2],Omega_M=0.3,Omega_Lambda=0.7,H0=1,\
             is_sn=1,err=1e-2,radi=5,grid=3,number_output=False):
        """
        Parameters
        ----------
        Omega_M : float , 
            component of visible mass. The default is 0.3.
        Omega_Lambda : float , 
            component of dark matter. The default is 0.7.
        Omega_K : float, 
            . The default is 0.0.
        H0 : float, optional
            present Hubble constant. The default is 1.
        ith : array,
            
        beta: array
            the source position according to lens

        Returns
        -------
        redshift of difference.
        """
        
        cosmo = cosmology.LambdaCDM(H0,Omega_M,Omega_Lambda)
        Dl = cosmo.angular_diameter_distance(zl)
        Ds = cosmo.angular_diameter_distance(zs)
        Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
        
        if is_sn:
            is_sn = 1e9
        else:
            is_sn = 1
        
        """
        value to ISs
        """
        
        
        #def ellipsoidal_lens(zl,zs,beta,ith=[]):
            
        images = self.solve_images(Dls,Ds,beta,err=err,radi=radi,grid=grid,\
                                   number_output=number_output)
        #print("number of images:",len(images))
        #print("images:",images)
        
        b = (4*Dls/Ds*self.sigma**2).decompose().value
        if type(b) != np.ndarray:
            b = np.array([b])
        
        if len(ith) == 0: #if ith is empty, gives out the whole redshift
            z_list = []
            for i in range(len(images)):
                
                z_temp = - (1+zl)*cosmo.H(zs)*\
                    (Dl*Ds)/(Dls*constants.c)\
                            *self.potential_tau(images[i],beta,to_arc=1)
                z_temp = z_temp.decompose().value
                z_list.append(z_temp)
                            
            return np.array(z_list)*is_sn
        
        elif len(ith) == 1: #if ith length is 1, gives ith redshift
            z_ith = - (1+zl)*cosmo.H(zs)*\
                (Dl*Ds)/(Dls*constants.c)\
                         *self.potential_tau(images[ith[0]],beta,to_arc=1)
            z_ith = z_ith.decompose().value
                         
            return z_ith*is_sn
                         
        elif len(ith) == 2: # if ith-lenth = 2, 
            Delta_z = - (1+zl)*cosmo.H(zs)*\
                    (Dl*Ds)/(Dls*constants.c)\
                        *(self.potential_tau(images[ith[0]],beta,to_arc=1) \
                             -self.potential_tau(images[ith[1]],beta,to_arc=1))
            Delta_z = Delta_z.decompose().value
            
            return Delta_z*is_sn
                            
        elif len(ith) == 3:
            z_list = []
            
            if len(b) == 1:
                for i in range(len(images)):
                    z_temp = - (1+zl)*cosmo.H(zs)*\
                        (Dl*Ds)/(Dls*constants.c)\
                                *self.potential_tau(images[i],beta,to_arc=1)
                    z_temp = z_temp.decompose().value
                    z_list.append(z_temp)
                    
                z_list = np.array([z_list])
            else:
                for i in range(len(images)):
                    z_li = []
                    self.b = b[i]
                    for j in range(len(images[i])):
                        z_temp = - (1+zl[i])*cosmo.H(zs[i])*\
                            (Dl[i]*Ds[i])/(Dls[i]*constants.c)\
                                    *self.potential_tau(images[i][j],beta,to_arc=1)
                        z_temp = z_temp.decompose().value
                        z_li.append(z_temp)
                    z_li = np.array(z_li)
                    z_list.append(z_li)
                    
                z_list = np.array(z_list)

            
            if len(z_list) == 1:
                Delta_z_max = z_list.max() - z_list.min()
                return np.array([Delta_z_max])
            
            elif len(z_list) >= 2:
                Delta_z_max = []
                for i in range(len(z_list)):
                    temp = z_list[i].max() - z_list[i].min()
                    Delta_z_max.append(temp)
            
                return np.array(Delta_z_max) #is_sn = 1e9
        
            
            
            
            
        
    def wCDM(self,zl,zs,beta,ith=[0,1,2],Omega_M=0.3,Omega_Lambda=0.7,\
             w0=-1,H0=1,is_sn=1,err=1e-2,radi=5,grid=3,number_output=False):
        """
        Parameters
        ----------
        Omega_M : float , 
            component of visible mass. The default is 0.3.
        Omega_Lambda : float , 
            component of dark matter. The default is 0.7.
        Omega_K : float, 
            . The default is 0.0.
        w0 : float
            the exponent of SoE, default value is -1
        H0 : float, optional
            present Hubble constant. The default is 1.
        ith : array,
            
        beta: array
            the source position according to lens 

        Returns
        -------
        redshift of difference.
        """
        
        cosmo = cosmology.wCDM(H0,Omega_M,Omega_Lambda,w0)
        Dl = cosmo.angular_diameter_distance(zl)
        Ds = cosmo.angular_diameter_distance(zs)
        Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
        
        if is_sn:
            is_sn = 1e9
        else:
            is_sn = 1
        
        """
        value to ISs
        """
        
        
        #def ellipsoidal_lens(zl,zs,beta,ith=[]):
            
        images = self.solve_images(Dls,Ds,beta,err=err,radi=radi,grid=grid,\
                                   number_output=number_output)
        #print("number of images:",len(images))
        #print("images:",images)
        
        b = (4*Dls/Ds*self.sigma**2).decompose().value
        if type(b) != np.ndarray:
            b = np.array([b])
        
        if len(ith) == 0: #if ith is empty, gives out the whole redshift
            z_list = []
            for i in range(len(images)):
                
                z_temp = - (1+zl)*cosmo.H(zs)*\
                    (Dl*Ds)/(Dls*constants.c)\
                            *self.potential_tau(images[i],beta,to_arc=1)
                z_temp = z_temp.decompose().value
                z_list.append(z_temp)
                            
            return np.array(z_list)*is_sn
        
        elif len(ith) == 1: #if ith length is 1, gives ith redshift
            z_ith = - (1+zl)*cosmo.H(zs)*\
                (Dl*Ds)/(Dls*constants.c)\
                         *self.potential_tau(images[ith[0]],beta,to_arc=1)
            z_ith = z_ith.decompose().value
                         
            return z_ith*is_sn
                         
        elif len(ith) == 2: # if ith-lenth = 2, 
            Delta_z = - (1+zl)*cosmo.H(zs)*\
                    (Dl*Ds)/(Dls*constants.c)\
                        *(self.potential_tau(images[ith[0]],beta,to_arc=1) \
                             -self.potential_tau(images[ith[1]],beta,to_arc=1))
            Delta_z = Delta_z.decompose().value
            
            return Delta_z*is_sn
                            
        elif len(ith) == 3:
            z_list = []
            
            if len(b) == 1:
                for i in range(len(images)):
                    z_temp = - (1+zl)*cosmo.H(zs)*\
                        (Dl*Ds)/(Dls*constants.c)\
                                *self.potential_tau(images[i],beta,to_arc=1)
                    z_temp = z_temp.decompose().value
                    z_list.append(z_temp)
                    
                z_list = np.array([z_list])
            else:
                for i in range(len(images)):
                    z_li = []
                    self.b = b[i]
                    for j in range(len(images[i])):
                        z_temp = - (1+zl[i])*cosmo.H(zs[i])*\
                            (Dl[i]*Ds[i])/(Dls[i]*constants.c)\
                                    *self.potential_tau(images[i][j],beta,to_arc=1)
                        z_temp = z_temp.decompose().value
                        z_li.append(z_temp)
                    z_li = np.array(z_li)
                    z_list.append(z_li)
                    
                z_list = np.array(z_list)

            
            if len(z_list) == 1:
                Delta_z_max = z_list.max() - z_list.min()
                return np.array([Delta_z_max])
            
            elif len(z_list) >= 2:
                Delta_z_max = []
                for i in range(len(z_list)):
                    temp = z_list[i].max() - z_list[i].min()
                    Delta_z_max.append(temp)
            
                return np.array(Delta_z_max) #is_sn = 1e9
            

    def w0waCDM(self,zl,zs,beta,ith=[0,1,2],Omega_M=0.3,Omega_Lambda=0.7,\
             w0=-1,wa=0,H0=1,is_sn=1,err=1e-2,radi=5,grid=3,number_output=False):
        """
        Parameters
        ----------
        Omega_M : float , 
            component of visible mass. The default is 0.3.
        Omega_Lambda : float , 
            component of dark matter. The default is 0.7.
        Omega_K : float, 
            . The default is 0.0.
        w0 : float
            the exponent of SoE, default value is -1
        H0 : float, optional
            present Hubble constant. The default is 1.
        ith : array,
            
        beta: array
            the source position according to lens 

        Returns
        -------
        redshift of difference.
        """
        
        cosmo = cosmology.w0waCDM(H0, Omega_M, Omega_Lambda,w0=w0,wa=wa)
        Dl = cosmo.angular_diameter_distance(zl)
        Ds = cosmo.angular_diameter_distance(zs)
        Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
        
        if is_sn:
            is_sn = 1e9
        else:
            is_sn = 1
        
        """
        value to ISs
        """
        
        
        #def ellipsoidal_lens(zl,zs,beta,ith=[]):
            
        images = self.solve_images(Dls,Ds,beta,err=err,radi=radi,grid=grid,\
                                   number_output=number_output)
        #print("number of images:",len(images))
        #print("images:",images)
        
        b = (4*Dls/Ds*self.sigma**2).decompose().value
        if type(b) != np.ndarray:
            b = np.array([b])
        
        if len(ith) == 0: #if ith is empty, gives out the whole redshift
            z_list = []
            for i in range(len(images)):
                
                z_temp = - (1+zl)*cosmo.H(zs)*\
                    (Dl*Ds)/(Dls*constants.c)\
                            *self.potential_tau(images[i],beta,to_arc=1)
                z_temp = z_temp.decompose().value
                z_list.append(z_temp)
                            
            return np.array(z_list)*is_sn
        
        elif len(ith) == 1: #if ith length is 1, gives ith redshift
            z_ith = - (1+zl)*cosmo.H(zs)*\
                (Dl*Ds)/(Dls*constants.c)\
                         *self.potential_tau(images[ith[0]],beta,to_arc=1)
            z_ith = z_ith.decompose().value
                         
            return z_ith*is_sn
                         
        elif len(ith) == 2: # if ith-lenth = 2, 
            Delta_z = - (1+zl)*cosmo.H(zs)*\
                    (Dl*Ds)/(Dls*constants.c)\
                        *(self.potential_tau(images[ith[0]],beta,to_arc=1) \
                             -self.potential_tau(images[ith[1]],beta,to_arc=1))
            Delta_z = Delta_z.decompose().value
            
            return Delta_z*is_sn
                            
        elif len(ith) == 3:
            z_list = []
            
            if len(b) == 1:
                for i in range(len(images)):
                    z_temp = - (1+zl)*cosmo.H(zs)*\
                        (Dl*Ds)/(Dls*constants.c)\
                                *self.potential_tau(images[i],beta,to_arc=1)
                    z_temp = z_temp.decompose().value
                    z_list.append(z_temp)
                    
                z_list = np.array([z_list])
            else:
                for i in range(len(images)):
                    z_li = []
                    self.b = b[i]
                    for j in range(len(images[i])):
                        z_temp = - (1+zl[i])*cosmo.H(zs[i])*\
                            (Dl[i]*Ds[i])/(Dls[i]*constants.c)\
                                    *self.potential_tau(images[i][j],beta,to_arc=1)
                        z_temp = z_temp.decompose().value
                        z_li.append(z_temp)
                    z_li = np.array(z_li)
                    z_list.append(z_li)
                    
                z_list = np.array(z_list)

            
            if len(z_list) == 1:
                Delta_z_max = z_list.max() - z_list.min()
                return np.array([Delta_z_max])
            
            elif len(z_list) >= 2:
                Delta_z_max = []
                for i in range(len(z_list)):
                    temp = z_list[i].max() - z_list[i].min()
                    Delta_z_max.append(temp)
            
                return np.array(Delta_z_max) #is_sn = 1e9
            
            
            
