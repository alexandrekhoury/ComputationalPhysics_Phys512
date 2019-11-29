# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 20:19:16 2019

@author: Alexandre
"""

import numpy as np
from matplotlib import pyplot as plt
plt.ion()

#this function is used for the error analysis and is presented in question2
def Ax(V,mask):
    #Vuse=np.zeros([V.shape[0]+2,V.shape[1]+2])
    #Vuse[1:-1,1:-1]=V
    Vuse=V.copy()
    Vuse[mask]=0
    ans=(Vuse[1:-1,:-2]+Vuse[1:-1,2:]+Vuse[2:,1:-1]+Vuse[:-2,1:-1])/4.0
    ans=ans-V[1:-1,1:-1]
    return ans
#size of grid box
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


V=bc.copy()

#calculating the real potential
"""
we know that at the edge of the box, the potential is 0 and 
that in the circle the potential is 1 (or at the edge of the circle)

the radius of the box is just half the size of a side of the box

we can calculate the difference in potential and find the constant that is
multiplying the expression (lambda). We can then isolate for the constant 
to allow ourselves to have a potential of 0 on the edges of the box and a
potential of 1 in the circle.
"""
V_box=0
V_circle=1
rbox=n//2
rcirc=rad

#finding the multiplier
V_tot=V_box-V_circle
lamb=V_tot/(np.log(rbox/rcirc))

#finding the constant (C is the same in both cases)
C=V[rbox+rcirc,rbox]-lamb*np.log(rcirc)
C=V[0,0]-lamb*np.log(rbox)

V_true=0.5*np.log((xv-nx//2)**2 +(yv-ny//2)**2)*lamb+C
V_true[mask]=bc[mask]
rho_true=V_true[1:-1,1:-1]-(V_true[1:-1,0:-2]+V_true[1:-1,2:]+V_true[:-2,1:-1]+V_true[2:,1:-1])/4.0


b=-(bc[1:-1,0:-2]+bc[1:-1,2:]+bc[:-2,1:-1]+bc[2:,1:-1])/4.0
tol=0.01
#calculating the numerical values 
for i in range(n**2):
    V[1:-1,1:-1]=(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0
    V[mask]=bc[mask]
#    plt.clf()
#    plt.imshow(V)
#    plt.colorbar()
#    plt.pause(0.001)
    r=b-Ax(V,mask)
    rtr=np.sum(r*r)
    
    if rtr<=tol:
        break    
print('on iteration ' + repr(i) + ' residual is ' + repr(rtr))    

rho=V[1:-1,1:-1]-(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0


#PLOTTING THE NUMERICAL POTENTIAL AND THE NUMERICAL DENSITY

plt.clf()
plt.imshow(V)
plt.colorbar()
plt.savefig("potential1")

plt.clf()
plt.imshow(rho)
plt.colorbar()
plt.savefig("density1")


#PLOTTING THE ANALYTICAL POTENTIAL AND THE ANALYTICAL DENSITY


plt.clf()
plt.imshow(V_true)
plt.colorbar()
plt.savefig("potential_true1")

plt.clf()
plt.imshow(rho_true)
plt.colorbar()
plt.savefig("density_true1")
    



