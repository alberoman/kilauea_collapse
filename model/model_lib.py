#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:35:16 2020

@author: aroman
"""

from scipy import optimize
import numpy as np

def Onechamber(V,beta,taustaud,Phi,a,l,mu,S,alpha,rhog):
    k = 1. / beta
    R5Samp =0               #Parameter inherited from old scaling of the equations. Must be equal to zero with the current scaling.    
    #Scaling
    Psi = rhog * V /(k*S)
    xstar = taustaud * V / (k * S)
    tstar =  8 * mu * l / (3.14 * a**4)
    dt = 0.05 #dt appproximately 20th of the shortest timescale
    N_cycles =int(((1 + Psi)/ (4 * alpha * Psi) * Phi))
    epsilon = 4 * alpha / Phi
    p_base = []
    x_base = []
    t_base = []
    tcounter = 0
    p0 =   4 * alpha /(1 + Psi)
    for i in range(1, N_cycles ):
        tslip =  np.log((N_cycles * (1 +  epsilon) - i) / (N_cycles - i ))
        N = int(tslip / dt)
        tsegment = np.linspace(0,tslip,N)
        psegment = (p0 + Phi) *  np.exp(-tsegment) - Phi
        x_base.append(4 * alpha / (1 + Psi) * i * np.ones(len((tsegment))))
        p_base.append(psegment)
        t_base.append(tsegment + tcounter)
        p0 =   + 4 * alpha / (1 + Psi) -4 * alpha * Psi /(1 + Psi) * i
        tcounter = tcounter + tslip
    x = np.concatenate(x_base)
    p = np.concatenate((p_base)) 
    t = np.concatenate((t_base)) 
    
   
    xPiston = x * xstar
    p = p * taustaud
    t = t * tstar
    ptpd = Phi * taustaud
    q = 3.14 * a**4 / (8 * mu * l) * (p + ptpd)
    Vcum = 0
    dt = np.diff(t) 
    VErupted = np.zeros(len(dt))
    for i in range(len(dt)):
        Vcum = Vcum + 0.5 * (q[i] + q[i+1]) * dt[i]
        VErupted[i] = Vcum
    
    return t,p,xPiston,VErupted
    

def p_feeder(t,R3,Gamma,Theta,a,b,c,d,pCollapse0,ps0):
    pFeeder = -R3 + (-R3 - pCollapse0 + (R3 + ps0)*(Gamma + a - Theta + 1)/2) * np.exp(t*(b - c))/a + (2*R3 + 2*a*(R3 + ps0) + 2*pCollapse0 - (R3 + ps0)*(Gamma + a - Theta + 1)) * np.exp(t*(b + c))/(2*a)
    return  pFeeder

def p_nonfeeder(t,R3,Gamma,Theta,a,b,c,d,pCollapse0,ps0):
    pNonFeeder = -R3 + (-c + d)*(-R3 - pCollapse0 + (R3 + ps0)*(Gamma + a - Theta + 1)/2) * np.exp(t*(b - c))/a + (c + d)*(2*R3 + 2*a*(R3 + ps0) + 2*pCollapse0 - (R3 + ps0)*(Gamma + a - Theta + 1)) * np.exp(t*(b + c))/(2*a)
    return pNonFeeder

def p_nonfeeder_root(t,R3,Gamma,Theta,a,b,c,d,pCollapse0,ps0,pslip):
    pNonFeederRoot = -R3 + (-c + d)*(-R3 - pCollapse0 + (R3 + ps0)*(Gamma + a - Theta + 1)/2) * np.exp(t*(b - c))/a + (c + d)*(2*R3 + 2*a*(R3 + ps0) + 2*pCollapse0 - (R3 + ps0)*(Gamma + a - Theta + 1)) * np.exp(t*(b + c))/(2*a) - pslip
    return  pNonFeederRoot


def segments_LF(w0,par,pslip,tslip,dt):
    pd0,pCollapse0 = w0
    r1,r3,r5,Gamma,Theta,a,b,c,d = par
    tslip = optimize.brentq(p_nonfeeder_root,0,1e+13, args = (r3,Gamma,Theta,a,b,c,d,pCollapse0,pd0,pslip))
    N = int(tslip / dt)
    tsegment = np.linspace(0,tslip,N)
    pdsegment = p_feeder(tsegment,r3,Gamma,Theta,a,b,c,d,pCollapse0,pd0)
    pssegment = p_nonfeeder(tsegment,r3,Gamma,Theta,a,b,c,d,pCollapse0,pd0)
    pdend = pdsegment[-1]
    return tsegment,tslip,pssegment,pdsegment,pdend


    
 

def Twochambers_LF(VDeep,VShallow,betaDeep,betaShallow,taustaud,PD0,Phi,ae,ac,le,lc,mu,S,alpha,rhog):
    kShallow = 1. / betaShallow
    kDeep = 1. / betaDeep
    R5Samp =0               #Parameter inherited from old scaling of the equations. Must be equal to zero with the current scaling.    
    #Scaling
    Psi = rhog * VShallow /(kShallow*S)
    Psi_bulk = rhog * (VShallow + VDeep) /(kShallow*S)
    Gamma = (ae /ac  )**4 * lc /le
    Theta = kShallow /kDeep * VDeep / VShallow
    tstar = VDeep * 8 * mu * lc / (kDeep * 3.14 * ac**4)
    xstar = taustaud * VShallow / (kShallow * S)
    A = np.sqrt(Gamma**2 - 2*Gamma*Theta + 2*Gamma + Theta**2 + 2*Theta + 1)
    B = -Gamma/2 - Theta/2 - 1./2
    C =  np.sqrt(4*Theta + (-Gamma + Theta - 1)**2)/2    
    D = Gamma/2 - Theta /2 + 1./2
    params = [Psi,Phi,R5Samp,Gamma,Theta,A,B,C,D]
    timescale1 = 1. / (A - B)
    timescale2 = 1. / (A + B)
    dt = min(timescale1,timescale2) / 20 #dt appproximately 20th of the shortest timescale
    pCollapse0 =   4 * alpha /(1 + Psi)
    PSLIP = - 4 * alpha * Psi * (1 - R5Samp)/(1 + Psi)
    TSLIP = 0
    w0 = np.array([PD0,pCollapse0])
    N_cycles =int(((1 + Psi)/ (4 * alpha * Psi) * Phi)) 
    i  = 1
    t_base = []
    ps_base = []
    pd_base = []
    x_base = []
    tcounter = 0
    for i in range(1,N_cycles + 1):
        tseg,tslip,PS,PD,PD0 = segments_LF(w0,params,PSLIP,TSLIP,dt)
        x_base.append(4 * alpha / (1 + Psi) * i * np.ones(len((tseg))))
        ps_base.append(PS)
        pd_base.append(PD)
        t_base.append(tseg + tcounter)
        pCollapse0 =   + 4 * alpha / (1 + Psi) -4 * alpha * Psi * (1 - R5Samp)/(1 + Psi) * i
        PSLIP =  - 4 * alpha * Psi * (1 - R5Samp)/(1 + Psi) * (i + 1)
        w0 = np.array([PD0,pCollapse0])
        tcounter = tcounter + tslip
    x = np.concatenate(x_base)
    pd = np.concatenate((pd_base)) 
    ps = np.concatenate((ps_base)) 
    t = np.concatenate((t_base)) 
    
   
    x = x * xstar
    ps = ps * taustaud
    pd = pd * taustaud
    t = t * tstar
    ptpd = Phi * taustaud
    q = 3.14 * ae**4 / (8 * mu * le) * (pd + ptpd)
    Vcum = 0
    dt = np.diff(t) 
    V = np.zeros(len(dt))
    
    for i in range(len(dt)):
        Vcum = Vcum + 0.5 * (q[i] + q[i+1]) * dt[i]
        V[i] = Vcum
    x = x[:-1]
    ps = ps[:-1]
    pd = pd[:-1]
    t = t[:-1]
    return t,ps,pd,x,V



