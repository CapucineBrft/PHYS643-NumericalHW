#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 11:39:42 2020

Advection Equation using the FTCS and the Lax-Friedrich methods.
@author: capucinebarfety

"""
import numpy as np
import matplotlib.pyplot as plt


#Define grid and time steps, as well as the input parameters of the functions
N = 250 #time increments
J = 50 #number of grid elements
dx = 1 
dt = 1

x = np.arange(1, J+1, step = dx) #spatial grid

u = -0.1 #velocity


#Create initial conditions
f_ftc = x.copy()*1./J
f_lf = x.copy()*1./J


#Set up the figure
plt.ion()
fig, axes = plt.subplots(1, 2, figsize=(10,6))


axes[0].plot(x, x/J, 'k-')
FTCS, = axes[0].plot(x,  f_ftc, linestyle='', marker='*')
axes[0].set_title('FTCS method')
axes[0].set_ylabel('f', fontsize=15)
axes[0].set_xlabel('x', fontsize=15)

axes[1].plot(x, x/J, 'k-')
LF, = axes[1].plot(x,  f_lf, linestyle='', marker='d')
axes[1].set_title('Lax-Friedrich method')
axes[1].set_ylabel('f', fontsize=15)
axes[1].set_xlabel('x', fontsize=15)

fig.canvas.draw()


#Iterate over time
for n in range(N):
    
    #for the FTC method:
    f_ftc[1:J-1] = f_ftc[1:J-1] - (u*dt/(2*dx)) * (f_ftc[2:] - f_ftc[:J-2])
    FTCS.set_ydata(f_ftc)
    
    
    #For the Lax-Friedrich method:
    f_lf[1:J-1] = 1./2. * (f_lf[2:] + f_lf[:J-2]) - (u*dt/(2*dx)) * (f_lf[2:] - f_lf[:J-2])
    LF.set_ydata(f_lf) 
    
    fig.canvas.draw()
    plt.pause(0.001)



    
    
    
    
    
    