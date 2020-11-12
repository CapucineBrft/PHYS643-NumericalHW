# PHYS643 Numerical Assignment
Capucine Barfety
Python 3.7

## Advection Problem:
The python script for this problem is "advection.py". 

The equation is solved using two methods, the Forward-Time Central-Space method and the Lax-Friedrich method. Looking at the figures created by the script, the FTCS solution is clearly unstable, as opposed to the Lax-Friedrich which is stable under the Courant Condition.



## Advection-Diffusion Problem:
The python script for this problem is "advection_diffusion.py".

The equation is solved using the implicit method for diffusion part, and the Lax-Friedrich method for the advection. Part of the code from the previous problem is used for the Lax-Friedrich method. The code defines $\delta t$ so that the Courant condition is respected, and the result is plotted for two different diffusion coefficients: D=0.01, and D=1.

As the velocity is constant, the boundary conditions for this problem are no-slip boundaries at both the right and left boundary of the grid, to ensure that there is no diffusive flux throught the first element.

We can see on the plots that the lower diffusion coefficient curve increases much faster than the other one. 


## 1D Hydro-Solver:
The python script for this problem is "HydroSolver.py".

The equations for density and momentum are solved using the donor cell advection method. 

The initial conditions for this problem are an initial constant density -  here set to 1 - to which we add a small Gaussian pertubation with ampltiude 0.01. Similarly, the velocity at N=0 (initial time) is defined as $v_0 = 0$ plus a Gaussian perturbation. As N increases, the perturbation appears to split into two gaussians that move in opposite directions towards both boundaries of the grid. There ot gets reflected back into a Gaussian.

When the amplitude of the perturbations is increased (to 0.3 for example), we can observe a shock when the to Gaussians meet in the original place: they form a spike that has a higher amplitude than the initial condition, and that is steeper. The width of this spike is set by the velocity. When we increase the amplitude of the perturbation in u, the spike gets sharper.
 
