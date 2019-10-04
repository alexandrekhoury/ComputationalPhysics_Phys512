# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 23:36:08 2019

@author: Alexandre
"""

import numpy as np
from matplotlib import pyplot as plt

time,flux,random =np.transpose(np.loadtxt('229614158_PDCSAP_SC6.txt',delimiter=','))

to=time[np.where(flux==np.max(flux))]#initial time
#isolating the region of interest
new_time=np.array(time[np.where(time==to)[0][0]:np.size(time)-500])
new_flux=np.array(flux[np.where(time==to)[0][0]:np.size(time)-500])

def func(params,t):
    #A is the amplitude
    #a is the rate of decay
    #to is the time-position of the flare
    #c is the initial y-axis value
    
    A=params[0]
    a=params[1]
    to=params[2]
    c=params[3]
    
    fun=A*np.exp(-a*(t-to))+c
    return t,fun

A=float(flux[np.where(flux==np.max(flux))]-1) #amplitude of the exponential
a=55
to=float(time[np.where(flux==np.max(flux))])#initial time
c=1
#guess parameters
guess1=np.array([A,a,to,c])
#the second and fourth value of the above array were find by trial and error


def grad_exp(params,t):
    
    
    A=params[0]
    a=params[1]
    to=params[2]
    c=params[3] 
       
    y=A*np.exp(-a*(t-to))+c
    
    grad=np.zeros([t.size,params.size])
    #now differentiate w.r.t. all the parameters
    
    grad[:,0]=np.exp(-a*(t-to))
    
    grad[:,1]=(-(t-to)*A*np.exp(-a*(t-to)))
    
    grad[:,2]=a*A*np.exp(-a*(t-to))
    
    grad[:,3]=t*0+1
    
    
    
    return y,grad

def newton(flux):
    for j in range(5):
        
    
        pred,grad=grad_exp(guess1,new_time)
        
        
        r=flux-pred
        err=(r**2).sum()
        r=np.matrix(r).transpose()
        grad=np.matrix(grad)
    
        lhs=grad.transpose()*grad
        rhs=grad.transpose()*r
        dp=np.linalg.inv(lhs)*(rhs)
        
        
        for jj in range(np.size(guess1)):
            guess1[jj]=guess1[jj]+float(dp[jj])
            
    return pred,guess1

def get_error():
    
    #adding random noise with rms amplitude then refitting using newton to see the uncertainty on each parameter
    
    pred,guess1=newton(new_flux)
    #finding the rms value to add random noise with rms amplitude
    Vrms=np.sqrt(np.mean((new_flux-pred)**2))
    
    new_guess=[]
    new_guess0=[]
    new_guess1=[]
    new_guess2=[]
    new_guess3=[]
    for i in range(7):
        random=np.random.normal(0,Vrms)
        noise=pred+random
        new_guess=(newton(noise)[1])
        
        new_guess0=np.append(new_guess0,new_guess[0])
        new_guess1=np.append(new_guess1,new_guess[1])
        new_guess2=np.append(new_guess2,new_guess[2])
        new_guess3=np.append(new_guess3,new_guess[3])

    err0=np.std(new_guess0)
    err1=np.std(new_guess1)
    err2=np.std(new_guess2)
    err3=np.std(new_guess3)
    
    err_array=([err0,err1,err2,err3])
    #these are the errors we get on each parameter
    return err_array
    
pred,guess1=newton(new_flux)

err=get_error()
print("The fit parameters are: "+str(guess1)+ "and the error is : " +str(err) )


    
plt.clf()  
plt.xlabel("Time")
plt.ylabel("Flux") 
plt.plot(new_time,new_flux);
plt.plot(new_time,func(guess1,new_time)[1]);
plt.plot(new_time,pred)
plt.legend(['Initial flux','Guess fit','Reduced best fit'])

