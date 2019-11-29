# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:22:44 2019

@author: Alexandre
"""

import numpy as np
from matplotlib import pyplot as plt
plt.ion()


def upres_mat(mat):
    """A quick & dirty way to increase the resolutio of a potential matrix by a factor
    of 2.  A smarter version here would lead to faster convergence, but I have ignored that."""
    mm=np.zeros([mat.shape[0]*2,mat.shape[1]*2],dtype=mat.dtype)
    mm[::2,::2]=mat
    mm[::2,1::2]=mat
    mm[1::2,::2]=mat
    mm[1::2,1::2]=mat
    return mm


def Ax(V,mask):

    Vuse=V.copy()
    Vuse[mask]=0
    ans=(Vuse[1:-1,:-2]+Vuse[1:-1,2:]+Vuse[2:,1:-1]+Vuse[:-2,1:-1])/4.0
    ans=ans-V[1:-1,1:-1]
    return ans

def pad(A):
    AA=np.zeros([A.shape[0]+2,A.shape[1]+2])
    AA[1:-1,1:-1]=A
    return AA

#initial resolution
n=128

V=np.zeros([n,n])  
rad=15

#increase the resolution by a facotr of 2, 6 times.
for i in range(6):
    if i >0:
        
        #increasing array size and radius size
        n=n*2
        rad=rad*2
    
    

    #making the cylinder
    nx, ny = (n, n)
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    xv, yv = np.meshgrid(x, y)
    
    
    
    #initializing the potential and the boundayry conditions
      
    bc=0*V
    
    
    #assigning the radius
    
    mask= (xv-nx//2)**2 +(yv-ny//2)**2 <=rad**2
    
    #getting the boundary conditions
    bc[mask]=1.0
    
    #masking the boundary of the masks
    mask[:,0]=True
    mask[:,-1]=True
    mask[0,:]=True
    mask[-1,:]=True
    
    
    
    b=-(bc[1:-1,0:-2]+bc[1:-1,2:]+bc[:-2,1:-1]+bc[2:,1:-1])/4.0
    
    #V=0*bc
    
    r=b-Ax(V,mask)
    p=r.copy()
    
    tol=0.01
    
    k=0
    
    #repeating the conjugate gradient method
    for k in range(n):
        Ap=(Ax(pad(p),mask))
        rtr=np.sum(r*r)
        
        alpha=rtr/np.sum(Ap*p)
    
        V=V+pad(alpha*p)
        rnew=r-alpha*Ap
        beta=np.sum(rnew*rnew)/rtr
        p=rnew+beta*p
        r=rnew
        
        
        
        
        
        
        
        if rtr<=tol:
            break
    print('on iteration ' + repr(k) + ' residual is ' + repr(rtr))
    #increasing the resolution
    
    plt.clf();
    plt.imshow(V)
    plt.colorbar()
    plt.savefig("res_"+str(n)+" X "+str(n) + "_potential3")
    if i!=5:
        V=upres_mat(V)
    
rho=V[1:-1,1:-1]-(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0


       

