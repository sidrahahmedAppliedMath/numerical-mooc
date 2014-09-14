# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 23:35:36 2014

@author: ken
"""

from math import sin, cos, log, ceil, pi
import numpy
import matplotlib.pyplot as plt
#%matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def f(u):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """
    
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    return numpy.array([-g*sin(theta) - C_D/C_L*g/v_t**2*v**2,
                      -g*cos(theta)/v + g/v_t**2*v,
                      v*cos(theta),
                      v*sin(theta)])
                      
def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    
    return u + dt * f(u)


def time_loop(u0,dt):
    u = numpy.zeros((1, 4))
    u[0] = u0 # fill 1st element with initial values
    # time loop - Euler method
    T = 0                               # final time
    
    y=u0[3]
    n=0
    while y >0:
    
        
        u = numpy.append(u,[euler_step(u[n], f, dt)],axis = 0)
        y = u[n+1,3]
        n=n+1
        T+=dt
    
    return u,T
    
# model parameters for paperplane:
g = 9.8      # gravity in m s^{-2}
v_t = 4.9   # trim velocity in m s^{-1}   
C_D = 1/5.  # drag coefficient --- or D/L if C_L=1
C_L = 1.0    # for convenience, use C_L = 1

### set initial conditions ###

x0 = 0.0  # horizotal position is arbitrary
y0 = 1.5  # initial altitude

v0_max = 10.5 #maximum initial velocity
thetha0_max = pi/3 #maximal initial angle

dv0 = 0.5 #velocity step
dt0 = 0.035 #angle step in radians

N_v0 = int(v0_max/dv0) + 1 
N_t0 = int(thetha0_max/dt0 ) + 1 

max_distances = numpy.zeros((N_v0,N_t0))
#ini_values = numpy.array([[i,j] for i in numpy.linspace(0.5,v0_max,N_v0) 
#for j in numpy.linspace(0,thetha0_max,N_t0)])
#
#
#triples = numpy.empty_like(ini_values,dtype=numpy.ndarray)
#for i,v0_theta0 in enumerate(ini_values):
v0_values = numpy.linspace(0.5,v0_max,N_v0) 
theta0_values = numpy.linspace(0,thetha0_max,N_t0)
for i in range(N_v0):
    v0 = v0_values[i]
    for j in range(N_t0):
        theta0 = theta0_values[j]
        
        
        
    

   
    
        dt = 0.001    
        u,T = time_loop(numpy.array([v0, theta0, x0, y0]),dt)
        max_distances[i,j]=  u[:,2][-1]


max_dist =  max_distances.max()    
max_ind = numpy.where( max_distances == max_dist )
max_ind = zip(max_ind[0],max_ind[1])[0]
v0 = v0_values[max_ind[0]]

theta0 = theta0_values[max_ind[1]] 

print "Maximum distance is %.2f. \
Initial values: v0 = %.2f, t0 = %.2f" % (max_dist,v0,theta0)



dt = 0.001    
u,T = time_loop(numpy.array([v0, theta0, x0, y0]),dt)
# get the glider's position with respect to the time
x = u[:,2]
y = u[:,3]
v = u[:,0]

# visualization of the path
plt.figure(figsize=(8,6))
plt.grid(True)
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'y', fontsize=18)
plt.title('Glider trajectory, flight time = %.2f , distance = %.2f' % (T,max(x)), fontsize=18)
plt.plot(x,y, 'k-', lw=2);

plt.figure(figsize=(8,6))
plt.pcolor(max_distances)
plt.colorbar()
plt.ylabel(r'$v_0$', fontsize=18)
plt.xlabel(r'$\theta_0$', fontsize=18)
plt.suptitle('Maximum distance (in m) for different initial conditions.', fontsize=18)
plt.show()
