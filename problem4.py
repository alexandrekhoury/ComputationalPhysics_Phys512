# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 20:53:46 2019

@author: Alexandre

All constants where set to 1 
"""

import numpy as np
from scipy import integrate

#defining the integrand of the electric field

def function(u,R,z):
    f=z-R*u
    f/=(R**2+z**2-2*R*z*u)**(3/2)

    return f

#defining the integral in scipy

def scipy_int(R,a,b):
    
    #making sure that R is in the values 
    z=np.linspace(R-5,R+4,1000)
    val=[]
    
    #taking values for all possible Z
    for i in z:
        
        #defining the function with u as the independant variable
        func= lambda u: function(u,R,i)
        
        #adding to the array a set of numerical results for different z of the integral
        val=np.append(val,integrate.quad(func,a,b)[0])
    plt.clf()    
    plt.plot(z,val)   
    return val


#taking the function in question3
def real_integrate(fun,first,middle,last,a,b,tol):
    
    x=np.linspace(a,b,5)
    
    dx=(b-a)/4.0
    
    
    #np.median(np.diff(x))
    
    #recalling the 1rst, 3rd and 5th point that are already evaluated (no need to do it twice)
    y=[first,fun(x[1]),middle,fun(x[3]),last]
    neval=2 #let's keep track of function evaluations
    f1=(y[0]+4*y[2]+y[4])/6.0*(b-a)
    f2=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12.0*(b-a)
    myerr=np.abs(f2-f1)
    #print([a,b,f1,f2])
    
    
    if (myerr<tol):
        #return (f2)/1.0,myerr,neval
        return (16.0*f2-f1)/15.0,myerr,neval
    else:
        mid=0.5*(b+a)
        f_left,err_left,neval_left=real_integrate(fun,y[0],y[1],y[2],a,mid,tol/2.0)
        f_right,err_right,neval_right=real_integrate(fun,y[2],y[3],y[4],mid,b,tol/2.0)
        neval=neval+neval_left+neval_right
        f=f_left+f_right
        err=err_left+err_right
        return f,err,neval


def integrate_Q3(R,a,b):
    

    
    #making sure that R is in the values 
    z=np.linspace(R-5,R+4,1000)
    val=[]
    
    #taking values for all possible Z
    for i in z:
        
        #defining the function with u as the independant variable
        func= lambda u: function(u,R,i)
        
        #if the try except error was not added, the integral would not work
        try:
        #adding to the array a set of numerical results for different z of the integral
            val=np.append(val,real_integrate(func,func(a),func((a+b)/2),func(b),a,b,1e-4)[0])
        except:
            val=np.append(val,0)
    plt.clf()    
    plt.plot(z,val)   
    return val

integrate_Q3(12,-1,1)


print("for the integrating method found in question 3, the singularity affects the function and could not resolve")

print("There seems to be a singularity at R. Scipy completely ignores it. The integrator does not care.")