# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 13:51:40 2019

@author: Alexandre
"""
import numpy as np
from matplotlib import pyplot as plt

#timestep
dt = 0.0005
#position step
dy = 0.0005
#constant of heat equation
k = 10**(-4)

#max position (box of length 1)
y_max = 0.04
#max time
t_max = 3
#initial temperature
T0 = 0

def FTCS(dt,dy,t_max,y_max,k,T0):
    
    #solving the discretized heat equation (putting all constants to the
    #right hand side of the equation) 
    s = k*dt/dy**2
    y = np.arange(0,y_max+dy,dy) 
    t = np.arange(0,t_max+dt,dt)
    r = len(t)
    c = len(y)
    T = np.zeros([r,c])
    
    #linear increase in temperature of the left side of the box
    T[:,0] = t
    for n in range(0,r-1):
        for j in range(1,c-1):
            
            #heat equation
            T[n+1,j] = T[n,j] + s*(T[n,j-1] - 2*T[n,j] + T[n,j+1]) 
        
        #non fixed boundary condition
        j = c-1 
        T[n+1, j] = T[n,j] + s*(T[n,j-1] - 2*T[n,j] + T[n,j-1])
    return y,T,r,s

y,T,r,s = FTCS(dt,dy,t_max,y_max,k,T0)

plot_times = np.arange(0.01,t_max,0.01)
for t in plot_times:

    plt.plot(y,T[int(t//dt),:])
    plt.ylim([0,T[:,].max()])
    
    #plt.pause(0.0001)
plt.savefig("heat_equation")   