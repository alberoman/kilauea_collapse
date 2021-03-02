#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 06:57:19 2021
@author: Alberto Roman (Jet Propulsion Laboratory)

This script calculates the main output of the caldera collapse model. It considers that the collapse proceeds in a shallow chamber (corresponding to the HLM source of Kilauea).
The main output are the pessure evolution in chamber, the total erupted volume, and the piston position as function of time.
Reservoirs pressure  are expressed relative to the critical shear stress to initiate motion (tau_d in the paper), so that if the pressure is 0 in the code output, it is actually equal to -tau_s. 
The pressure timeseries can be used to calculate surface deformation with any deformation source model. One should be careful to consider that the bulk compressibilties (beta parameters) are
beta = beta_m + beta_c
where beta is the bulk compressibilty beta_m is the compressibilty of the magma, and beta_c is the compessibility of the reservoir (which depends on its geometrical features). One should specify thus the bulk compressibility accordingly with the derformation source model chosen.
"""

from model_lib import *

#Fixed parameters
rho = 2700      #Density
g = 9.8         #Gravity  
rhog = rho * g


#Conduit parameters
l = 4e+4               #Length of the conduit feeding the eruption (meters)
a = 3                  #Width of the conduit feeding the eruption (meters)
#Pressures and friction parameters
taustaud = 1e+6 #Difference between static and dynamic friction (Pa)
phi = 20    #Dimensionless parameters topography pressure to friction stress (dimensionless). Note that this parameters should be always greater than one to generate several collapses. Range of this parameter  5-100 seems reasonable. Typically the higher the parameter the larger the number of collapse events. 

#Magma viscosity
mu = 5 # (Pa s)

#Cylinder parameters
Rcyl = 5e+2   #Radius (m)
S = 3.14 * Rcyl**2 
alpha = 1   #Aspect ratio (Dimensionless)

#Reservoirs volume (V) and bulk compressibilities (beta). S
beta = 1e-9 # (Pa**-1)
V = 1e+10   #(m**3)

t,p,xPiston,VErupted = Onechamber(V,beta,taustaud,phi,a,l,mu,S,alpha,rhog)  
#output is time expressed in seconds (t), the pressure in the reservoir expressed in pascals,piston position (in meters) the volume of the erupted material (in m**3)  

