# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 19:09:58 2019

@author: Alexandre
"""

import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0,2,1000)

y1=np.exp(x)
y2=np.exp(0.01*x)



#for e^x
delta1= 10**(-16/5)

#for e^(0.01x)
delta2=(10**-16/(0.01)**(5))**(1/5)

print("The delta value for e^x is :" + str(delta1))
print("this is of order 10^-3")
print("")
print("The delta value for e^x is :" + str(delta2))
print("this is of order 10^-1")