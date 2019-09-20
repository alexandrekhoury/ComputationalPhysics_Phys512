# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 12:08:00 2019

@author: Alexandre

I interpolated for values of V, then projected thoses values on the x-axis to get the temperature.
This should logically still give similar results.

All 
"""

import numpy as np
from scipy.interpolate import splrep,splev
import matplotlib.pyplot as plt

a=0.09
b=1.65
steps=500

#loading data 
T,V,D =np.transpose(np.loadtxt("lakeshore.txt"))

#starting the interpolation
tck = splrep(V[::-1] ,T[::-1])

Vnew=np.linspace(a,b,steps)

Tnew=splev(Vnew,tck)
#end of interpolation


#usefull for error calculation (having same number of elements in matrix)
Terr=splev(V,tck)


#now if we want to get a value, we would call in from the command line a description 
#as follows: get_temp(V) where V is any voltage
    
def get_temp(V):
    
    #gets the temperature from a given Voltage

    temp=Tnew[np.where(V<=Vnew)[-1][-1]]

    return temp

def get_error_rough():
    
    #I am taking the interpolated values and substracting the real values,
    #then i take the mean to get an average error
    #gets the rough estimate of the error which is found to be of order 10^-13
    err=np.max(Terr-T)
    print(str(err))


get_error_rough()

def get_error_less_rough(f,f_true):
    #taking the derivatives at each point and comparing them 
    #using the formula in question 1
    
    f_der=[]
    f_der_true=[]
    for i in range (2,len(f)-2):
        f_der=np.append(f_der,(8*(f[i+1]-f[i-1])-f[i+2]+f[i-2])/(12*(V[i+1]-V[i])))
        f_der_true=np.append(f_der_true,(8*(f_true[i+1]-f_true[i-1])-f_true[i+2]+f_true[i-2])/(12*(V[i+1]-V[i])))
    
    print(np.max(f_der_true-f_der))

get_error_less_rough(Terr,T)

#calling the above function gets the rough estimate of the error in the derivatives (order 10^-15) in the first terms
#taking the mean we get an error once again of order 10^-12    

def even_better_error(Ti):
    #interpolate again with less data points
    a=0.09
    b=1.65
    steps=500
    
    
    #taking out a random data point in the interpolation, then evaluating the error there
    
    T_err=np.delete(T[::-1],34)
    V_err=np.delete(V[::-1],34)
    tck = splrep(V_err,T_err )

    Vnew_err=np.linspace(a,b,steps)

    Tnew_err=splev(Vnew_err,tck)


    print(np.max(Ti-Tnew_err))
    

even_better_error(Tnew)  
#the error is of order 10^-6 which is already a better estimate but still very small
#in the above function, i randomly removed a point and interpolated again
  


plt.clf()
plt.plot(Vnew,Tnew)
plt.plot(V,T,marker=".")