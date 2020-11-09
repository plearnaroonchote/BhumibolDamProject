#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 22:56:44 2018

@author: plearnaroonchote
"""

import numpy as np
import math

def Simulate(storage_t, release, inflow, evaporation, MAX_S, unit, spill,
             pump, MAX_R):
    storage_t1 = storage_t + inflow*unit - release*unit-evaporation*unit
    + pump*unit - spill*unit;
    if (storage_t1 < 0):
        release = (0-storage_t1)*(1/unit);
        spill = 0
        return 0, release, spill;
    elif (storage_t1 > MAX_S):
        release = MAX_R;
        spill = (storage_t1 - MAX_S)*(1/unit) - MAX_R;
        return MAX_S, release, spill;
    else: 
        return storage_t1, release, spill;
    
def Power_Generation(power_release,efficiency,forebay_elev,tailwater_elev):
    rel_rate = power_release*(10**6)/(24*60*60)
    if(rel_rate > 1200):
        P= 1200*9.81*(forebay_elev-tailwater_elev)*efficiency
    else:
        P= rel_rate*9.81*(forebay_elev-tailwater_elev)*efficiency
    return P;  
 
def Objectives(release, demand, spill, power):   
    obj1 = np.max(release + spill);
    obj2 = np.sum(np.square(np.maximum(0,(demand - release))));
    obj3 = np.sum(power);
    return obj1, obj2, obj3;

def Elevation(storage):
    ForebayElevation = (-3*10**-7)*(storage**2) + 0.0098*storage + 180.51 
    return ForebayElevation

def LinearOperatingCurve(S, Hf, D):   
    if S > Hf:
        R = D
    else:
        m = D/Hf
        R = m*(S)
    return R

def NewCubicRadialBasis(vector, weight, center, radius, num_vars):
    
    vector[-2] = np.sin(2*math.pi*vector[:,-1]/365) #how do we calculate the phase shift??
    vector[-1] = np.cos(2*math.pi*vector[:,-1]/365)
    
    value = np.empty(num_vars)
    
    for var in range(num_vars):
        for i in range(3):
           weight[var,i]= weight[var,i]/np.sum(weight[var,:])
           value[var] +=  weight[var,i] * (abs((vector[var] - center[var,i])/radius[var,i])**3)
        
    return value;
    
def CubicRadialBasis(vector, weight, center, radius):
    
    value = 0;
    
    for i in range(3):
          weight[i]= weight[i]/np.sum(weight)
          value +=  weight[i] * (abs((vector - center[i])/radius[i])**3)
        
    return value;
    
def Forecasts(Qt,Xt, Q, rho_qf):
    
    Rqx = 0.94290158
    X_mean = np.mean(Q[:-1])
    Q_mean = np.mean(Q[1:])
    F_mean = np.mean(Q)
    Sq = np.std(Q[1:])
    Sx = np.std(Q[:-1])
    Sigf = Sq * rho_qf 
    Sigf2 = (rho_qf**2) * (Sq**2)
    d2 = (Sigf2/(Sx**2)) * ((1-(rho_qf**2))/(1-(Rqx**2)))
    
    d = np.sqrt(d2)
    b = rho_qf*(Sigf/Sq) - d*Rqx*(Sx/Sq)
    F = F_mean + b*(Qt-Q_mean) + d*(Xt-X_mean)
    
    return F;
    
    
    

    
    
