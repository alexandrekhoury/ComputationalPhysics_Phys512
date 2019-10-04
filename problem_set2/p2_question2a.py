# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 21:14:05 2019

@author: Alexandre
"""

import numpy as np
import matplotlib.pyplot as plt



time,flux,random =np.transpose(np.loadtxt('229614158_PDCSAP_SC6.txt',delimiter=','))

to=time[np.where(flux==np.max(flux))]#initial time
#isolating the region of interest
new_time=time[np.where(time==to)[0][0]:np.size(time)]
new_flux=flux[np.where(time==to)[0][0]:np.size(time)]

def func(params,t):
    #A is the amplitude
    #a is the rate of decay
    #to is the time-position of the flare
    #c is the initial y-axis value
    
    A=params[0]
    a=params[1]
    to=params[2]
    c=params[3]
    t=t[np.where(t==to)[0][0]:np.size(t)]
    fun=A*np.exp(-a*(t-to))+c
    return t,fun

A=flux[np.where(flux==np.max(flux))]-1 #amplitude of the exponential
a=55
to=time[np.where(flux==np.max(flux))]#initial time
c=1
#guess parameters
guess=[A,a,to,c]
#the second and fourth value of the above array were find by trial and error

plt.clf()
plt.xlabel("Time")
plt.ylabel("Flux")
plt.xlim(1706.4,1706.7)
plt.plot(time,flux)
x,y=func(guess,time)
plt.plot(x,y)
plt.savefig("Q2_a_exponential")

