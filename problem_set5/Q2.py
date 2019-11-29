# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 12:46:00 2019

@author: Alexandre
"""


# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 20:19:16 2019

@author: Alexandre
"""

import numpy as np
from matplotlib import pyplot as plt
plt.ion()


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

#grid size
n=400

#initializing arrays
V=np.zeros([n,n])
bc=0*V

#making the cylinder
nx, ny = (n, n)
x = np.linspace(0, nx-1, nx)
y = np.linspace(0, ny-1, ny)
xv, yv = np.meshgrid(x, y)
#assigning the radius
rad=50
mask= (xv-nx//2)**2 +(yv-ny//2)**2 <=rad**2

#getting the boundary conditions
bc[mask]=1.0

#masking the boundary of the masks
mask[:,0]=True
mask[:,-1]=True
mask[0,:]=True
mask[-1,:]=True


#setting the right hand side of the equation Ax=b
b=-(bc[1:-1,0:-2]+bc[1:-1,2:]+bc[:-2,1:-1]+bc[2:,1:-1])/4.0

V=0*bc

r=b-Ax(V,mask)
p=r.copy()

#setting the tolerance
tol=0.01

#repeating the process untill the tolerance has been reached
for k in range(n):
    Ap=(Ax(pad(p),mask))
    rtr=np.sum(r*r)
    
    alpha=rtr/np.sum(Ap*p)

    V=V+pad(alpha*p)
    rnew=r-alpha*Ap
    beta=np.sum(rnew*rnew)/rtr
    p=rnew+beta*p
    r=rnew
    plt.clf();
    plt.imshow(V)
    plt.colorbar()
    plt.pause(0.0001)
    
    if rtr<=tol:
        break
print('on iteration ' + repr(k) + ' residual is ' + repr(rtr))    
rho=V[1:-1,1:-1]-(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0

plt.clf()
plt.imshow(V)
plt.colorbar()
plt.savefig("potential2")
