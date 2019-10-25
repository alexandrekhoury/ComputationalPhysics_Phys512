# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 23:03:21 2019

@author: Alexandre
"""

import numpy as np
from matplotlib import pyplot as plt
import camb

wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
datax,datay,err=[wmap[:,0],wmap[:,1],wmap[:,2]]



def get_chisquare(model,data,err):
    chisq=np.sum((model-data)**2/err**2)
    return chisq

def covar_mtx(func,err):
    noise = np.diag(err**2)
    cov = np.linalg.inv(np.dot(np.dot(func.transpose(),noise), func))
    perr = np.sqrt(np.diag(cov))
    
    return cov,perr

def get_spectrum(par,lmax):
    #print('pars are ',pars)
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

def take_step(scale_factor=1.0):
    param_scalings=np.asarray([1,0.01,0.05,0.01,1e-9,0.1])
    return scale_factor*param_scalings*np.random.randn(len(param_scalings))

def take_step_cov(covmat):
    mychol=np.linalg.cholesky(covmat)
    return np.dot(mychol,np.random.randn(covmat.shape[0]))


def MC_run(datay,err):
    params=np.asarray([65,0.02,0.1,0.05,2e-9,0.96])
    
    #params=pars+take_step()
    
    nstep=5000
    chains_new=np.zeros([nstep,len(params)]) #keep track of where the chain went
    size=datax.size
    cmb=get_spectrum(params,wmap[:,0].size)[2:]
    chisq=get_chisquare(cmb,datay,err)
    scale_fac=0.5
    #COVARIANCE MATRIX FROM QUESTION 2
    mycov=np.asarray([[ 1.24868057e+01,  1.98089582e-03, -2.20099477e-02,
         3.71676590e-01,  1.37176141e-09,  7.46404593e-02],
       [ 1.98089582e-03,  6.03960807e-07, -2.55027928e-06,
         8.19655275e-05,  3.24253016e-13,  1.74000853e-05],
       [-2.20099477e-02, -2.55027928e-06,  4.40462208e-05,
        -6.09360124e-04, -2.16622908e-12, -1.15536164e-04],
       [ 3.71676590e-01,  8.19655275e-05, -6.09360124e-04,
         2.09550017e-02,  8.14721363e-11,  3.00779694e-03],
       [ 1.37176141e-09,  3.24253016e-13, -2.16622908e-12,
         8.14721363e-11,  3.18596617e-19,  1.16578953e-11],
       [ 7.46404593e-02,  1.74000853e-05, -1.15536164e-04,
         3.00779694e-03,  1.16578953e-11,  6.09937697e-04]])
    
    new_params=params+take_step_cov(mycov)*scale_fac
    chisqvec_new=np.zeros(nstep)
    for i in range(nstep):
        
        new_params=params+take_step_cov(mycov)*scale_fac
        if (new_params[3]>0):
            
            new_model=get_spectrum(new_params,size)[2:]
            new_chisq=get_chisquare(new_model,datay,err)
        
            delta_chisq=new_chisq-chisq
            prob=np.exp(-0.5*delta_chisq)
            accept=np.random.rand(1)<prob
            if accept:
                params=new_params
                chisq=new_chisq 
            chains_new[i,:]=params
            chisqvec_new[i]=chisq
        
    fit_params=np.mean(chains_new,axis=0)
    
    return chains_new,fit_params,chisqvec_new


chains,fit,chisqvec=MC_run(datay,err)


chains2=chains.copy()

#getting rid of zeroes terms
for i in range (0,np.where(chains2[:,0]<10)[0].size):
    chains2=np.delete(chains2,np.where(chains2[:,0]<10)[0][0],0)

fit2=np.zeros(chains2[0,:].size)
std=np.zeros(chains2[0,:].size)
for i in range (0,chains2[0,:].size):
    fit2[i]=np.mean(chains2[int(nstep/2):,i])
    std[i]=np.std(chains2[int(nstep/2):,i])
    
for i in range (1,7):
    plt.subplot(610+i)
    plt.plot(chains2[:,i-1])
    
"""   
#for fft plots 
for i in range (1,7):
    plt.subplot(610+i)
    plt.plot(np.fft.rfft(chains2[:,i-1])[10:chains2[:,i-1].size])
"""  
#to find the error : 
mycov=np.cov(chains2[int(chains2[:,1].size/2):,:].T)
error=np.sqrt(np.diag(mycov))

print("The error on the parameters are:"+ str(error))




plt.ion()