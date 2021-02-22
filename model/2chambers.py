#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 12:17:22 2020

@author: aroman
"""

import matplotlib.pyplot as plt
import numpy as np
from model_lib import *
import numpy as np
import pickle

path_results = '../../../results/'

#Fixed parameters
rho = 2700
g = 9.8
rhog = rho * g


#Conduit parameters
le = 4e+4
lc = 2.5e+3
ae = 4 
ac = 4
#Sources volumes and compressibilities

#Pressures and friction6
pspd = 1e+6
phi = 10
#Magma viscosity
mu = 100


#Cylinder parameters
Rcyl = 1e+3
S = 3.14 * Rcyl**2 
#Reservoirs volume and bulk compressibilities
beta = 1e-9
Vdeep = 1e+10
Vcollapse= 1e+9

tv,volume,t,ps,pd,q,tcrit = Twochambers_syntethic_LF(Vdeep,Vcollapse,beta,pspd,phi,ae,ac,le,lc,mu,S,rhog)
