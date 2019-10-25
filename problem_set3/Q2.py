# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:09:25 2019

@author: Alexandre
"""

import numpy as np
import camb
from matplotlib import pyplot as plt
import time

wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
datax,datay,err=[wmap[:,0],wmap[:,1],wmap[:,2]]
def covar_mtx(func,err):
    noise = np.diag(1/err**2)
    cov = np.linalg.inv(np.dot(np.dot(func.transpose(),noise), func))
    perr = np.sqrt(np.diag(cov))
    
    return cov,perr

#for everyparameter add the small change dx in order to take the derivative
def change(pars,i,add):
    p=pars.copy()
    p[i]=pars[i]+add
    
    return p

def get_derivative(pars,i,size):
    #taking the derivatives at each point and comparing them 
    #using the formula in question 1 in assignment 1

    dx=pars[i]/10000
    
    #long derivative since takes more time to evaluate each function
    #f_der=(8*(get_spectrum(change(pars,i,dx),size)[2:]-get_spectrum(change(pars,i,-dx),size)[2:]-get_spectrum(change(pars,i,2*dx),size)[2:]+get_spectrum(change(pars,i,-2*dx),size)[2:]))/(12*(dx))
    #shorter derivative
    f_der=(get_spectrum(change(pars,i,dx),size)[2:]-get_spectrum(change(pars,i,-dx),size)[2:])/(2*(dx))
    
    return (f_der)

#use this function for fixed tau
def get_spectrum(par,lmax):
    #print('pars are ',pars)
    
    H0=par[0]
    ombh2=par[1]
    omch2=par[2]
    
    tau=0.05
    As=par[3]
    ns=par[4]
    
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    pars.InitPower.set_params(As=As,ns=ns,r=0)
    pars.set_for_lmax(lmax,lens_potential_accuracy=0)
    results=camb.get_results(pars)
    powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
    cmb=powers['total']
    tt=cmb[:,0]    #you could return the full power spectrum here if you wanted to do say EE
    
    return tt

def get_spectrum2(par,lmax):
    #print('pars are ',pars)
    #THIS FUNCTION IS ONLY USED TO GET THE COVARIANCE MATRIX FOR ALL 6 PARAMETERS AS MENTIONNNED IN QUESTION 3
    
    
    H0=par[0]
    ombh2=par[1]
    omch2=par[2]
    
    tau=par[3]
    As=par[4]
    ns=par[5]
    
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    pars.InitPower.set_params(As=As,ns=ns,r=0)
    pars.set_for_lmax(lmax,lens_potential_accuracy=0)
    results=camb.get_results(pars)
    powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
    cmb=powers['total']
    tt=cmb[:,0]    #you could return the full power spectrum here if you wanted to do say EE
    
    return tt

def grad_exp(pars,size):
    

    grad=np.zeros([size,pars.size])
    #now differentiate w.r.t. all the parameters
    p=pars.copy()

    for i in range(0,p.size):   
        grad[:,i]=get_derivative(p,i,size)
           
    return grad

def newton(data,pars,size,error):
    
    p=pars.copy()
    for j in range(5):
    
        pred=get_spectrum(p,size)[2:]
        grad=grad_exp(p,size)
        cov=covar_mtx(grad,error)
        r=data-pred
        err=(r**2).sum()
        r=np.matrix(r).transpose()
        grad=np.matrix(grad)
    
        lhs=grad.transpose()*grad
        rhs=grad.transpose()*r
        dp=np.linalg.inv(lhs)*(rhs)
        
        
        for jj in range(np.size(p)):
            p[jj]=p[jj]+dp[jj]

    
    return p,grad,cov

plt.ion()


pars=np.asarray([65,0.02,0.1,2e-9,0.96]) #fixed tau
#pars=np.asarray([65,0.02,0.1,0.05,2e-9,0.96]) #non fixed tau


#function that returns the chisquare given the model, data and variance
def get_chisquare(model,data,err):
    chisq=np.sum((model-data)**2/err**2)
    return chisq

cmb=get_spectrum(pars,wmap[:,0].size)[2:]
params,grad,cov=newton(datay,pars,datax.size,err)

#get the Chisquare for the new parameters
cmb2=get_spectrum(params,wmap[:,0].size)[2:]

chisquare=get_chisquare(cmb2,datay,err)
print("The error on the parameters are:" + str(cov[1]))

plt.clf();
#plt.errorbar(wmap[:,0],wmap[:,1],wmap[:,2],fmt='*')
#plt.plot(wmap[:,0],wmap[:,1],'.')


#plt.plot(cmb)

#FOR PLOTTING

plt.subplot(1, 2, 1)
plt.plot(cmb)
plt.title('CMB')



plt.subplot(1, 2, 2)
plt.plot(grad[:,0])
plt.title("Derivative for first parameter")

#calling the above function gets the rough estimate of the error in the derivatives (order 10^-15) in the first terms
#taking the mean we get an error once again of order 10^-12    