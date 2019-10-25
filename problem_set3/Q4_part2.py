# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 12:05:58 2019

@author: Alexandre
"""

import numpy as np
import math


chains=np.loadtxt("chains_Q3.txt")

#weighing the tau parameters according to a gaussian distribution
def gauss(x, x0, sigma):
    return 1/(2*np.pi)**(1/2)*1/sigma*np.exp(-(x-x0)**2/(2*sigma**2))

tau=0.0544
sigma=0.0073


weight=np.zeros(chains[:,3].size)
for i in range(0,chains[:,3].size):
    weight[i]=gauss(chains[:,3][i],tau,sigma)

params=np.zeros(chains[0,:].size)
std=np.zeros(chains[0,:].size)
for i in range(0,chains[0,:].size):
    params[i]=np.average(chains[:,i],weights=weight)
    std[i]=10*np.sqrt(np.average((chains[:,i]-(np.zeros((chains[:,i].size))+params[i]))**2, weights=weight))
    
    
    
print("The new parameters are:"+str(params))
print("The errors of the parameters " + str(std))