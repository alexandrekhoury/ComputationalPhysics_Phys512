# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:51:10 2019

@author: Alexandre
"""

import numpy as np
import matplotlib.pyplot as plt





def func(x):
    
    #defining a function for log base 2
    fun=np.log2(x)
    
    return fun

def cheb_mtx(nx,ord):
    
    #defining the chebychev polynomials for a certain order
    
    x=np.linspace(-1,1,nx)
    
    mat=np.zeros([nx,ord+1])
    
    mat[:,0]=1.0
    
    #assigning the chebychev polynomials to the matrix
    if ord>0:
        mat[:,1]=x
    
    if ord>1:
        for i in range(1,ord):
        
            mat[:,i+1]=2*x*mat[:,i]-mat[:,i-1]
            
    return mat,x

def apply_cheb():
    
    n=5000
    ord=150
    ncoeff=7
    mat,x=cheb_mtx(n,ord)
    
    x=np.linspace(0.5,1,n)
    
    #chebyshev
    lhs=np.dot(mat.transpose(),mat)
    rhs=np.dot(mat.transpose(),func(x))
    fitp=np.dot(np.linalg.inv(lhs),rhs)
    pred=np.dot(mat[:,:ncoeff],fitp[:ncoeff])
    
    
    #least squares:
    lhs2=np.dot(mat[:,:ncoeff].transpose(),mat[:,:ncoeff])
    rhs2=np.dot(mat[:,:ncoeff].transpose(),func(x))
    fitp2=np.dot(np.linalg.inv(lhs2),rhs2)
    pred2=np.dot(mat[:,:ncoeff],fitp2)
    

    plt.clf()
    plt.xlabel('X')
    plt.ylabel("Residuals")
    plt.plot(x,pred-func(x));plt.plot(x,pred2-func(x));

    plt.legend(['Chebyshev Residual','Least-Squares Residual'])
    
    plt.savefig("Q1_a_Chebyshev.png")

    
    print('rms error for chebyshev     is ',np.sqrt(np.mean((pred-func(x))**2)),' with max error ',np.max(np.abs(pred-func(x))))
    print('rms error for least-squares is ',np.sqrt(np.mean((pred2-func(x))**2)),' with max error ',np.max(np.abs(pred2-func(x))))
    print("7 coefficients are needed for the chebyshev polynomial to get an error less than or equal to 1e-6")
    

    
if __name__ == "__main__":
    apply_cheb()