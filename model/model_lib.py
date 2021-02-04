#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:35:16 2020

@author: aroman
"""

import numpy as np
from scipy import optimize
from scipy import interpolate
import matplotlib.pyplot as plt
import time 

def ps_analytic(t,R3,T1,phi,a,b,c,d,pd0,ps0):
    ps = -R3 + (-R3 - pd0 + (R3 + ps0)*(T1 + a - phi + 1)/2) * np.exp(t*(b - c))/a + (2*R3 + 2*a*(R3 + ps0) + 2*pd0 - (R3 + ps0)*(T1 + a - phi + 1)) * np.exp(t*(b + c))/(2*a)
    return  ps

def pd_analytic(t,R3,T1,phi,a,b,c,d,pd0,ps0):
    pd = -R3 + (-c + d)*(-R3 - pd0 + (R3 + ps0)*(T1 + a - phi + 1)/2) * np.exp(t*(b - c))/a + (c + d)*(2*R3 + 2*a*(R3 + ps0) + 2*pd0 - (R3 + ps0)*(T1 + a - phi + 1)) * np.exp(t*(b + c))/(2*a)
    return pd

def pd_analytic_root(t,R3,T1,phi,a,b,c,d,pd0,ps0,pslip):
    ps_root = -R3 + (-c + d)*(-R3 - pd0 + (R3 + ps0)*(T1 + a - phi + 1)/2) * np.exp(t*(b - c))/a + (c + d)*(2*R3 + 2*a*(R3 + ps0) + 2*pd0 - (R3 + ps0)*(T1 + a - phi + 1)) * np.exp(t*(b + c))/(2*a) - pslip
    return  ps_root


def TwoChambers_LF(w0,par,pslip,tslip):
    ps0,pd0 = w0
    r1,r3,r5,t1,phi,a,b,c,d = par
    tslip = optimize.brentq(pd_analytic_root,0,1e+13, args = (r3,t1,phi,a,b,c,d,pd0,ps0,pslip))
    tsegment = np.linspace(0,tslip,50)
    pssegment = ps_analytic(tsegment,r3,t1,phi,a,b,c,d,pd0,ps0)
    pdsegment = pd_analytic(tsegment,r3,t1,phi,a,b,c,d,pd0,ps0)
    psend = pssegment[-1]
    return tsegment,tslip,pssegment,pdsegment,psend


    
 

def Twochambers_syntethic_LF(VsSamp,VdSamp,betadSamp,pspdSamp,R3Samp,condsSamp,conddSamp,ls,ld,mu,S,rhog):
    alphaSamp = 1
    kdSamp = 1 / betadSamp 
    R5Samp =0
    ksSamp = kdSamp
    R1Samp = rhog * VdSamp /(kdSamp*S)
    T1 = (condsSamp / conddSamp )**4 * ld /ls
    PHI = kdSamp /kdSamp * VsSamp / VdSamp
    params = [T1,PHI,R3Samp] 
    tstar = VsSamp * 8 * mu * ld / (ksSamp * 3.14 * conddSamp**4)
    xstar = pspdSamp * VdSamp / (kdSamp * S)

    A = np.sqrt(T1**2 - 2*T1*PHI + 2*T1 + PHI**2 + 2*PHI + 1)
    B = -T1/2 - PHI/2 - 1./2
    C =  np.sqrt(4*PHI + (-T1 + PHI - 1)**2)/2    
    D = T1/2 - PHI /2 + 1./2
    params = [R1Samp,R3Samp,R5Samp,T1,PHI,A,B,C,D]
    
    PD0 =   4 * alphaSamp /(1 + R1Samp)
    PSLIP = - 4 * alphaSamp * R1Samp * (1 - R5Samp)/(1 + R1Samp)
    TSLIP = 0
    PS0 = 0
    w0 = np.array([PS0,PD0])
    TSLIP = 0
    N_cycles =int(((1 + R1Samp)/ (4 * alphaSamp * R1Samp) * R3Samp))
    i  = 1
    t_base = []
    ps_base = []
    pd_base = []
    tcounter = 0
    Ncrit = int(0.6 * N_cycles)
    print(N_cycles)
    for i in range(1,N_cycles + 1):
        tseg,tslip,PS,PD,ps0 = TwoChambers_LF(w0,params,PSLIP,TSLIP)
        ps_base.append(PS)
        pd_base.append(PD)
        t_base.append(tseg + tcounter)
        PD0 =   + 4 * alphaSamp / (1 + R1Samp) -4 * alphaSamp * R1Samp * (1 - R5Samp)/(1 + R1Samp) * i
        PS0 = ps0
        PSLIP =  - 4 * alphaSamp * R1Samp * (1 - R5Samp)/(1 + R1Samp) * (i + 1)
        w0 = np.array([PS0,PD0])
        tcounter = tcounter + tslip
        if i == Ncrit:
            tcrit = tcounter

    pd = np.concatenate((pd_base)) 
    ps = np.concatenate((ps_base)) 
    t = np.concatenate((t_base)) 

   
    
    ps = ps * pspdSamp
    pd = pd * pspdSamp
    t = t * tstar
    tcrit = tcrit * tstar
    ptpd = R3Samp * pspdSamp
    q = 3.14 * condsSamp**4 / (8 * mu * ls) * (ps + ptpd)
    Vcum = 0
    dt = np.diff(t) 
    V = np.zeros(len(dt))
    tVolume = np.zeros(len(dt))
    
    for i in range(len(dt)):
        Vcum = Vcum + 0.5 * (q[i] + q[i+1]) * dt[i]
        tVolume[i] = 0.5 * (t[i] + t[i + 1])
        V[i] = Vcum
    tVolume = tVolume
    return tVolume,V,t,ps,pd,q,tcrit



