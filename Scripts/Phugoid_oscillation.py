# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 19:16:43 2014

@author: Pierre
"""

import numpy
import matplotlib.pyplot as plt
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
for n in range(1,N):
    u = u + dt*numpy.array([u[1],g*(1-u[0]/zt)])
    z[n] = u[0]

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
plt.show()

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
    for n in range(1,N):
        ### compute next solution using Euler method ###
        u = u + dt*numpy.array([u[1], g*(1-u[0]/zt)])
        z[n] = u[0]   # store the elevation at time-step n+1
    
    z_values[i] = z.copy()    # store the total elevation calculation grid i

