#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 12:17:22 2020
@author: Alberto Roman (Jet Propulsion Laboratory)

This script calculates the main output of the two chamber model for caldera collapse. It considers that the collapse proceeds in a shallow chamber (corresponding to the HLM source of Kilauea), whereas the eruption is fed by the deeper chamber (corresponding to SCR for Kilauea).
The main output are the pessure evolution in the two chambers, the total erupted volume, and the piston position as function of time.
Reservoirs pressure  are expressed relative to the critical shear stress to initiate motion (tau_d in the paper), so that if the pressure is 0 in the code output, it is actually equal to -tau_s. 
The pressure timeseries can be used to calculate surface deformation with any deformation source model. One should be careful to consider that the bulk compressibilties (beta parameters) are
beta = beta_m + beta_c
where beta is the bulk compressibilty beta_m is the compressibilty of the magma, and beta_c is the compessibility of the reservoir (which depends on its geometrical features). One should specify thus the bulk compressibility accordingly with the derformation source model chosen.
The model allow two specify a compressibility for the upper reservoir (betaShallow) and one for the lower one (betaDeep) to account different geometries.
"""

from model_lib import *

#Fixed parameters
rho = 2700      #Density
g = 9.8         #Gravity  
rhog = rho * g


#Conduit parameters
le = 4e+4               #Length of the conduit feeding the eruption (meters)
lc = 2.5e+3             #Length of the conduit connecting the 2 chambers (meters)
ae = 4                  #Width of the conduit feeding the eruption (meters)
ac = 4                  #Width of the conduit cnnecting the 2 chambers (meters)

#Pressures and friction parameters
taustaud = 1e+6 #Difference between static and dynamic friction (Pa)
phi = 5    #Dimensionless parameters topography pressure to friction stress (dimensionless). Note that this parameters should be always greater than one to generate several collapses. Range of this parameter  5-100 seems reasonable. Typically the higher the parameter the larger the number of collapse events. 

#Magma viscosity
mu = 100 # (Pa s)

#Cylinder parameters
Rcyl = 5e+2   #Radius (m)
S = 3.14 * Rcyl**2 
alpha = 1   #Aspect ratio (Dimensionless)

#Reservoirs volume (V) and bulk compressibilities (beta). Shallow corresponds to the reservoir in which there is collapse (corresponds to subscript s in the paper). Deep corresponds to the reservoir feeding the eruption (corresponds to subscript s in the paper)
betaShallow = 1e-9 # (Pa**-1)
betaDeep = 1e-9    # (Pa**-1)   
VDeep = 1e+10   #(m**3)
VShallow= 1e+9  #(m**3)
pDeep0 = 0   #Initial pressure in the deep reservoir (Pa). 0 corresponds to pDeep = - taus. Pressures represent the difference to the lithostatic value. THe initial pressure in the shallow reservoirs does not have to be set, as it is equal to the pressure increase due to the firtst collapse
t,pShallow,pDeep,xPiston,VErupted = Twochambers_LF(VDeep,VShallow,betaDeep,
                                           betaShallow,taustaud,pDeep0,
                                           phi,ae,ac,le,lc,mu,S,alpha,rhog)  
#output is time expressed in seconds (t), the pressure in the shallow and deep reservoirs expressed in pascals,piston position (in meters) the volume of the erupted material (in m**3)  

