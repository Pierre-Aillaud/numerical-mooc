# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 19:16:43 2014

@author: Pierre
"""

import numpy
import matplotlib.pyplot as plt

def euler(u,z,dt):
    """Use Euler's method to solve the system of two first order ODE"""
    N = len(z)
    for n in range(1,N):
     u = u + dt*numpy.array([u[1],g*(1-u[0]/zt)])
     z[n] = u[0]
    return z
    
def get_error(z, dt):
    """Returns the error relative to analytical solution using L-1 norm.
    
    Parameters
    ----------
    z : array of float
        numerical solution.
    dt : float
        time increment.
        
    Returns
    -------
    err : float
        L_{1} norm of the error with respect to the exact solution.
    """
    N = len(z)
    t = numpy.linspace(0.0, T, N)
    
    z_exact = v*(zt/g)**.5*numpy.sin((g/zt)**.5*t)+\
                (z0-zt)*numpy.cos((g/zt)**.5*t)+zt
    
    return dt * numpy.sum(numpy.abs(z-z_exact))
    
    
#Grid in time
T = 100
dt = 0.01
N = int(T/dt)+1
t = numpy.linspace(0.0,T,N)

#Initial Condition
z0 = 100.
v = 10
g = 9.81
zt = 100.

u = numpy.array([z0, v])
z = numpy.zeros(N)
z[0] = z0

#time-loop using Euler's method
euler(u,z,dt)

#Exact solution
z_exact = v*(zt/g)**.5*numpy.sin((g/zt)**.5*t)+\
            (z0-zt)*numpy.cos((g/zt)**.5*t)+zt  
    
#plot the solution
plt.figure(figsize=(10,6))
plt.ylim(40,160)
plt.xlabel('t',fontsize=14)
plt.ylabel('z', fontsize=14)
plt.plot(t,z,'b-')
plt.plot(t,z_exact,'k-')
plt.legend(['Numerical solution', 'Analytical solution'],numpoints=1)

# time-increment array
dt_values = numpy.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0001])

# array that will contain solution of each grid
z_values = numpy.empty_like(dt_values, dtype=numpy.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T/dt)+1    # number of time-steps
    ### discretize the time using numpy.linspace() ###
    t = numpy.linspace(0.0, T, N)

    # initial conditions
    u = numpy.array([z0, v])
    z = numpy.empty_like(t)
    z[0] = z0
    
    # time loop - Euler method
    euler(u,z,dt)
    
    z_values[i] = z.copy()    # store the total elevation calculation grid i
    
error_values = numpy.empty_like(dt_values)


for i, dt in enumerate(dt_values):
    ### call the function get_error() ###
    error_values[i] = get_error(z_values[i], dt)
    
plt.figure(figsize=(10, 6))
plt.tick_params(axis='both', labelsize=14) #increase tick font size
plt.grid(True)                         #turn on grid lines
plt.xlabel('$\Delta t$', fontsize=16)  #x label
plt.ylabel('Error', fontsize=16)       #y label
plt.loglog(dt_values, error_values, 'ko-')  #log-log plot
plt.axis('equal')                      #make axes scale equally;
    
plt.show()

