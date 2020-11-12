#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:13:37 2020

Advection-Diffusion equation using the Lax-Friedrich and the implicit method.
@author: capucinebarfety

"""

import numpy as np
import matplotlib.pyplot as plt

###Using the same initial and boundary conditions as the advection problem:
###Define time step, as well as the input parameters of the functions
N = 75
u = -0.1
dx = 1
dt = 0.5*dx/np.abs(u)

#Diffusion coefficients:
D1 = 0.01
D2 = 1

beta1 = (D1*dt) /((dx)**2)
beta2 = (D2*dt) /((dx)**2)

#spatial grid:
x = np.arange(1, (N+1)*dx, step = dx)


###Define the tri-diagonal matrix
def matrix(b):
    a = (1+ 2*b)*np.eye(N) + np.eye(N, k=1)* -b  + np.eye(N, k=-1)*-b
    a[0][0]=1
    a[0][1]=0
    
    a[-1][-1]=1
    a[-1][-2]=0
    
    return a
    
A1 = matrix(beta1)
A2 = matrix(beta2)


###Initial boundary conditions
f1 = x.copy()/(N)
f2 = x.copy()/(N)


###Set up the figure
plt.ion()
fig, axes = plt.subplots(1,2, figsize=(15,6))
#Set plot for D1
axes[0].plot(x,f1, 'k-')
axes[0].set_title(f'D = {D1}')
axes[0].set_xlabel('x', fontsize=15)
axes[0].set_ylabel('f', fontsize=15)
AD1, = axes[0].plot(x, f1, linestyle='', marker='.')
#Set plot for D2
axes[1].plot(x,f2, 'k-')
axes[1].set_title(f'D = {D2}')
axes[1].set_xlabel('x', fontsize=15)
axes[1].set_ylabel('f', fontsize=15)
AD2, = axes[1].plot(x, f2, linestyle='', marker='.')

fig.canvas.draw()


###Iterate over time
for n in range(N):
    
    A1 = matrix(beta1)
    A2 = matrix(beta2)
    
    #Solving the equation for D1:
    f1 = np.linalg.solve(A1, f1)
    f1[1:N-1] = 1./2. * (f1[2:] + f1[:N-2]) - ( u*dt/(2*dx)) * (f1[2:] - f1[:N-2])
    AD1.set_ydata(f1) 
    
    #Solving the equation for D2:
    f2 = np.linalg.solve(A2, f2)
    f2[1:N-1] = 1./2. * (f2[2:] + f2[:N-2]) - ( u*dt/(2*dx)) * (f2[2:] - f2[:N-2])
    AD2.set_ydata(f2)
    
    fig.canvas.draw()
    plt.pause(0.0001)