#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 20:14:25 2020

1D Hydro-Solver using the donor cell advection solution
@author: capucinebarfety
"""
import numpy as np
import matplotlib.pyplot as plt

J = 100 #number of grid elements
N = 250 #time increments
x = np.arange(J)/(J-1) * (J) #spatial grid
dx = np.abs(x[1]-x[0])
dt = 0.15
cs = 5. #speed of sound

###Set the perturbations
#Create NxJ 0 matrices:
f1 = np.zeros((N, J)) 
f2 = np.zeros((N, J))
m = 0.5*(x[0]+x[-1]) #mean of the Gaussians
sigma = 0.1*J #sigma for the Gaussians
#Set the N=0 (so initial condition) of the matrices using density/velocity with gaussian perturbations
f1[0,:] = 1+0.01*np.exp( -(x-m)**2 / sigma**2) #rho_0 = 1
f2[0,:] = (1+0.01*np.exp( -(x-m)**2 / sigma**2))*(0.01*np.exp( -(x-m)**2 / sigma**2)) #v_0


###Set up the figure
plt.ion()
fig, axes = plt.subplots(1,1, figsize=(15,6))
axes.plot(x,f1[0], 'k-')
den, = axes.plot(x, f1[0], linestyle='', marker='d')
axes.set_xlabel('x', fontsize=15)
axes.set_ylabel(r'$\rho$', fontsize=15)
fig.canvas.draw()


###Iterate over time
for n in range(1,N):
    
    #Compute the velocity
    u = np.zeros(J+1)
    for i in range(1, J):
        u[i] = 0.5 *  (f2[n-1, i]/f1[n-1, i] + f2[n-1, i-1]/f1[n-1, i-1])
    
    #Compute J1    
    J1 = np.zeros(J+1)
    for i in range(1, J):
        if u[i] > 0:
            J1[i] = f1[n-1, i-1]*u[i]
        else:
            J1[i] = f1[n-1, i]*u[i]

    #Update f1:
    for i in range(J):
        f1[n, i]= f1[n-1, i] - dt/dx * (J1[i+1] - J1[i])


    #Compute J2:
    J2 = np.zeros(J+1)
    for i in range(1,J):
        if u[i] > 0:
            J2[i] = f2[n-1, i-1]*u[i]
        else:
            J2[i] = f2[n-1, i]*u[i]

    #Update f2:        
    for i in range(J):
        f2[n, i] = f2[n-1, i] - dt/dx * (J2[i+1] - J2[i])
    

    #Take into account pressure:
    for i in range(1, J-1):
        f2[n, i] = f2[n, i] - dt/dx * cs * (f1[n, i+1] - f1[n, i-1])
    
    
    #Reflective Boundary conditions
    f1[n, 0] = f1[n, 0] - (dt / dx) * J1[0] 
    f1[n, -1] = f1[n, -1] + (dt / dx) * J1[-2]

    f2[n, 0] = f2[n, 0] - (dt/ dx) * J2[0] 
    f2[n, -1] = f2[n, -1] + (dt / dx) * J2[-2]
    
    #plot the density
    den.set_ydata(f1[n])
    axes.set_title(f'Density at N={n}', fontsize=20)
    fig.canvas.draw()
    plt.pause(0.001)
    

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_ylabel('Time', fontsize=15)
ax.set_xlabel('x', fontsize=15)
ax.set_zlabel(r'$\rho$', fontsize=15)
ax.set_title('3D plot of the density as function of time and space', fontsize=20)
X,Y = np.meshgrid(x, np.arange(1,N+1))
ax.plot_surface(X, Y, Z = f1)
