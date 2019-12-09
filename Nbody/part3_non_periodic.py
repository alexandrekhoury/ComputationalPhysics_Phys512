# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:01:15 2019

@author: Alexandre
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 21:10:26 2019

@author: Alexandre
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class particles:
    def __init__(self,m=1.0,npart=10,soft=0.01,G=1.0,dt=0.1,grid_size=50):
        self.opts={}
        self.opts['soft']=soft
        self.opts['n']=npart
        self.opts['G']=G
        self.opts['dt']=dt
        self.opts['grid_size']=grid_size

        self.opts['m']=m
        self.x,self.y=np.random.randint(0, self.opts['grid_size'], size=(2, self.opts['n']))
        #self.m=m
        self.vx=np.double(self.x.tolist())*0
        self.vy=np.double(self.y.tolist())*0
        xx= np.linspace(0, grid_size-1, grid_size)


        kx,ky=np.meshgrid(xx,xx)
        
        soft=1
        Gr=1/(1e-13+4*np.pi*((kx)**2+(ky)**2)**(1/2))
        Gr[0,0]=1/(4*np.pi*soft)
        
        Gr+=np.flip(Gr,0)
        Gr+=np.flip(Gr,1)
        Gr_ft=np.fft.fft2(Gr)
        self.Gr_ft=Gr_ft

        self.energy=0
    def get_grid(self,x,y):
        
        grid_size=self.opts['grid_size']
    
        A=np.histogram2d((np.round(x)).astype(int),(np.round(y)).astype(int),bins=grid_size,range=[[0, grid_size], [0, grid_size]])[0]*self.opts['m']
        #A=np.zeros((self.opts['grid_size'],self.opts['grid_size']))
    
#        for i in range(0,self.opts['n']):
#            
#            
#            A[int(round(x[i]))%self.opts['grid_size'],int(round(y[i]))%self.opts['grid_size']]+=self.opts['m']         
#    
#    
        Gr_ft=self.Gr_ft
        rho_ft=np.fft.fft2(A)
    
        conv=Gr_ft*rho_ft
        pot=np.fft.ifft2(conv)
        pot=0.5*(np.roll(pot,1,axis=1)+pot)
        pot=0.5*(np.roll(pot,1,axis=0)+pot)
        
        return A,pot
        
    def get_force(self,x,y,pot,A):
        
        forcex=-1/2*(np.roll(pot,1,axis=0)-np.roll(pot,-1,axis=0))*A
        forcey=-1/2*(np.roll(pot,1,axis=1)-np.roll(pot,-1,axis=1))*A
        
        
        x_new=np.double(x.tolist())*0
        y_new=np.double(y.tolist())*0
        
        
        self.vx+=np.real(forcex[(np.round(x)%self.opts['grid_size']).astype(int),(np.round(y)%self.opts['grid_size']).astype(int)])*self.opts['dt']
        x_new=x+self.vx*self.opts['dt']
        self.vy+=np.real(forcey[(np.round(x)%self.opts['grid_size']).astype(int),(np.round(y)%self.opts['grid_size']).astype(int)])*self.opts['dt']
        y_new=y+self.vy*self.opts['dt']
        
        m=self.opts['m']
        self.energy=1/2*np.sum(self.vx**2+self.vy**2)*m+np.sum(pot)/2

        return x_new,y_new

        
if __name__=='__main__':
    

    dt=1
    n=100000
    grid_size=500
    part=particles(m=1,npart=n,dt=dt,grid_size=grid_size)
    time=500

    A,pot=part.get_grid(part.x,part.y)
    x_new,y_new=part.get_force(part.x,part.y,pot,A)
    
    grid=np.zeros([int(time//dt),grid_size,grid_size])
    count=0
    
    energy=np.zeros(int(time//dt))
    
    plt.figure()
    plt.pause(15)
    
    for i in np.arange(0,time,part.opts['dt']):
        
        
        A_new,pot=part.get_grid(x_new,y_new)
        x_new,y_new=part.get_force(x_new,y_new,pot,A_new)
        
#        if count<energy.shape[0]:
#            grid[count]=A_new
        energy[count]=np.real(part.energy)
        
        plt.clf()
        plt.imshow(abs(A_new))
        plt.pause(0.0001)
        
        count+=1
    

    fig, ax = plt.subplots(figsize=(8, 6))
    cax = ax.imshow(grid[0])
    
    cb = fig.colorbar(cax)
    
    def animate(i):
          
        cax.set_array(grid[i])
    
        
    anim = FuncAnimation(fig, animate, interval=10, frames=grid.shape[0], repeat=True,blit=False,save_count=grid.shape[0])
    
    anim.save("nparticles_nonperiodic.gif")