# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 14:13:05 2019

@author: Alexandre
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:22:44 2019

@author: Alexandre
"""

import numpy as np
from matplotlib import pyplot as plt


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
    #Vuse=np.zeros([V.shape[0]+2,V.shape[1]+2])
    #Vuse[1:-1,1:-1]=V
    Vuse=V.copy()
    Vuse[mask]=0
    ans=(Vuse[1:-1,:-2]+Vuse[1:-1,2:]+Vuse[2:,1:-1]+Vuse[:-2,1:-1])/4.0
    ans=ans-V[1:-1,1:-1]
    return ans

def pad(A):
    AA=np.zeros([A.shape[0]+2,A.shape[1]+2])
    AA[1:-1,1:-1]=A
    return AA


n=128
V=np.zeros([n,n])  
rad=30
z=5
for i in range(z):
    if i >0:
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
    #adding the bump
    mask2=(xv-nx//2)**2 +(yv-ny//2-rad)**2 <=(0.1*rad)**2
    mask[mask2]=True
    
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
    for k in range(n):
        Ap=(Ax(pad(p),mask))
        rtr=np.sum(r*r)
        
        alpha=rtr/np.sum(Ap*p)
    
        V=V+pad(alpha*p)
        rnew=r-alpha*Ap
        beta=np.sum(rnew*rnew)/rtr
        p=rnew+beta*p
        r=rnew
        
        
        
        
        
#        plt.clf();
#        plt.imshow(V)
#        plt.colorbar()
#        plt.pause(0.0001)
        
        if rtr<=tol:
            break
        
    print('on iteration ' + repr(k) + ' residual is ' + repr(rtr))
    if i!=z-1 :  
        V=upres_mat(V)
        
rho=V[1:-1,1:-1]-(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0 
plt.ion()

#different way of plotting electric field
#E=np.array(np.gradient(V))
#Elec=(E[0]**2+E[1]**2)**(0.5)
#
#
#fig, ax = plt.subplots() 
#plt.imshow(V)    
#ax.streamplot(x, y, -E[1], -E[0], linewidth=1, cmap=plt.cm.inferno,
#              density=2, arrowstyle='->', arrowsize=1.5)       
#        
#        

#reduce electric field computations so that the graph isn't all black by arrows
#only compute 1 electric field value for every 30 (size) potential values
size=30
V1=V[::size,::size]      
E=np.array(np.gradient(V1))
Elec=(E[0]**2+E[1]**2)**(0.5)


fig, ax = plt.subplots()

plt.imshow(V)
ax.quiver(xv[::size,::size],yv[::size,::size],-E[1],E[0])
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_aspect('equal')
plt.savefig("elec_fied_and_potential4")

#ELECTRIC FIELD AT BUMP:

E=np.array(np.gradient(V))
Elec=(E[0]**2+E[1]**2)**(0.5)

plt.clf()
plt.plot(Elec[n//2,:])
plt.savefig("elec_bump")