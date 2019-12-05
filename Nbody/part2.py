# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:30:07 2019

@author: Alexandre
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class particles:
    def __init__(self,m=1.0,npart=10,soft=1,G=1.0,dt=0.1,grid_size=50):
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
        
        Gr=1/(1e-13+4*np.pi*((kx)**2+(ky)**2)**(1/2))
        Gr[0,0]=1/(4*np.pi*soft)
        
        Gr+=np.flip(Gr,0)
        Gr+=np.flip(Gr,1)
        Gr_ft=np.fft.fft2(Gr)
        self.Gr_ft=Gr_ft
        
        self.energy=0

        
    def get_grid(self,x,y):
        
        grid_size=self.opts['grid_size']
    
        A=np.histogram2d((np.round(x)%grid_size).astype(int),(np.round(y)%grid_size).astype(int),bins=grid_size,range=[[0, grid_size], [0, grid_size]])[0]*self.opts['m']
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
        #making sure the potential is centered with the particle positions
        pot=0.5*(np.roll(pot,1,axis=1)+pot)
        pot=0.5*(np.roll(pot,1,axis=0)+pot)
        
        return A,pot
        
    def get_force(self,x,y,pot,A):
        
        #taking the gradient of the potential to get the forces
        #multiplying by density matrix
        
        forcex=-1/2*(np.roll(pot,1,axis=0)-np.roll(pot,-1,axis=0))*A
        forcey=-1/2*(np.roll(pot,1,axis=1)-np.roll(pot,-1,axis=1))*A
        
        #initializing new x and y positions 
        x_new=np.double(x.tolist())*0
        y_new=np.double(y.tolist())*0
        
        
        self.vx+=np.real(forcex[(np.round(x)%self.opts['grid_size']).astype(int),(np.round(y)%self.opts['grid_size']).astype(int)])*self.opts['dt']
        x_new=x+self.vx*self.opts['dt']
        self.vy+=np.real(forcey[(np.round(x)%self.opts['grid_size']).astype(int),(np.round(y)%self.opts['grid_size']).astype(int)])*self.opts['dt']
        y_new=y+self.vy*self.opts['dt']
        
        

        m=self.opts['m']
        self.energy=1/2*np.sum(self.vx**2+self.vy**2)*m+np.sum(pot)/2
        return x_new,y_new
    def evolve(self):
        
        kinetic=0.5*numpy.sum(self.m*(self.vx**2+self.vy**2))
        return pot+kinetic
        
if __name__=='__main__':
    
    
    

    n=2
    grid_size=150
    m=20
    dt=1
    part=particles(m=m,npart=n,dt=dt,grid_size=grid_size)
    time=405
    
    #initializing orbit (restricting initial position and velocity)
    part.x[0]=grid_size//2-15
    part.x[1]=grid_size//2+1
    
    part.y[0]=grid_size//2-15
    part.y[1]=grid_size//2+1
    
    part.vx[0]=0
    part.vx[1]=0
    
    part.vy[0]=1
    part.vy[1]=-1
    
    A,pot=part.get_grid(part.x,part.y)
    x_new,y_new=part.get_force(part.x,part.y,pot,A)
    
    
    """
    The commented code in the section below is to provide the videos.

    
    """
    

#    grid=np.zeros([int(time//dt),grid_size,grid_size])
#    count=0
    
    for i in np.arange(0,time,part.opts['dt']):
        
        
        A_new,pot=part.get_grid(x_new,y_new)
        
#        if count<grid.shape[0]:
#            grid[count]=A_new
        x_new,y_new=part.get_force(x_new,y_new,pot,A_new)
        
        plt.clf()
        plt.imshow(abs(A_new))
        plt.pause(0.0001)
#        count+=1
        

#    fig, ax = plt.subplots(figsize=(8, 6))
#    cax = ax.imshow(grid[0])
#    
#    cb = fig.colorbar(cax)
#    
#    def animate(i):
#          
#        cax.set_array(grid[i])
#    
#        
#    anim = FuncAnimation(fig, animate, interval=40, frames=4050, repeat=True,blit=False,save_count=4050)
#    
#    anim.save("2particles.gif")