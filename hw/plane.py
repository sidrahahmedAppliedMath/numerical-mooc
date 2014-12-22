# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 17:01:28 2014

@author: Ken
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 23:35:36 2014

@author: ken
"""

from math import sin, cos, log, ceil, pi, degrees
import numpy
import matplotlib.pyplot as plt
from scipy import optimize
#%matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'arial'
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

def euler_midpoint_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method with 
    middle point.
    
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
    u_star = u + 0.5*dt * f(u)
    return u + dt * f(u_star)

def time_loop(u0,dt):
    """Returns the array of u vectors at each time step
    
    Parameters
    ----------
    u0 : array of float
         initial conditions
    dt : float
         time step
         
    Returns
    -------
    u : array of u arrays with shape (n,4), n - number of time steps while y >0
    T : float
        final time
    """
    
    u = numpy.zeros((1, 4))
    u[0] = u0 # fill 1st element with initial values
    # time loop - Euler method
    T = 0                               # final time
    
    y=u0[3]
    n=0
    while y >0:
    
#        u = numpy.resize(u,(n+2,4))
#        u[n+1] = euler_step(u[n], f, dt)
        u = numpy.append(u,[euler_midpoint_step(u[n], f, dt)],axis = 0)
        y = u[n+1,3]
        n=n+1
        T+=dt
    
    return u,T

def time_loop_for_cycle(u0,dt):
    """Returns the u vector at final time step
    
    Parameters
    ----------
    u0 : array of float
         initial conditions
    dt : float
         time step
         
    Returns
    -------
    u : array of float
    T : float
        final time
    """
    u = numpy.copy(u0) # fill 1st element with initial values
    # time loop - Euler method
    T = 0                               # final time
    
    y=u0[3]
    n=0
    while y >0:
    
        u = euler_midpoint_step(u, f, dt)
#        u = numpy.append(u,[euler_step(u[n], f, dt)],axis = 0)
        y = u[3]
        n=n+1
        T+=dt
    
    return u,T
def distance(v0_theta0, *params):
    """Returns the distance that paperplane flies before hitting the ground
    when launched at given theta_0
    
    Parameters
    ----------
    u0 : array of float
         initial conditions
    
    Returns
    -------
    l : float
        final distance
    """
    v0,theta0 = v0_theta0
    x0,y0,dt = params
    
    u,T = time_loop_for_cycle(numpy.array([v0, theta0, x0, y0]),dt)
#    if u[2]<= 0:
#        print u[2]
#        dt = 0.1*dt
#        u,T = time_loop_for_cycle(numpy.array([v0, theta0, x0, y0]),dt)
#        print u[2]
    return u[2]
def inv_distance(v0_theta0, *params):
    v0,theta0 = v0_theta0
    x0,y0,dt = params
    return 1./distance(v0_theta0, *params)    
def neg_distance(v0_theta0, *params):
    v0,theta0 = v0_theta0
    x0,y0,dt = params
    return -distance(v0_theta0, *params)    


dt = 0.001
# model parameters for paperplane:
g = 9.8      # gravity in m s^{-2}
v_t = 4.9   # trim velocity in m s^{-1}   
C_D = 1/5.  # drag coefficient --- or D/L if C_L=1
C_L = 1.0    # for convenience, use C_L = 1

### set initial conditions ###

x0 = 0.0  # horizotal position is arbitrary
y0 = 2  # initial altitude

v0_max = 12.5 #maximum initial velocity
thetha0_max = pi/3 #maximal initial angle

dv0 = 0.5 #velocity step
dt0 = 0.035 #angle step in radians

N_v0 = int(v0_max/dv0) + 1 
N_t0 = int(thetha0_max/dt0 ) + 1 

max_distances = numpy.zeros((N_v0,N_t0))

v0_values = numpy.linspace(0.5,v0_max,N_v0) 
theta0_values = numpy.linspace(-0.5,thetha0_max,N_t0)

rranges = (slice(0.5,v0_max,dv0),slice(-0.5,thetha0_max,dt0))
resbrute = optimize.brute(neg_distance, rranges, args=(x0,y0,dt), full_output=True,
                          finish=optimize.fmin)

#resbrute_none = optimize.brute(neg_distance, rranges, args=(x0,y0,dt), full_output=True,
#                          finish=None)

#==============================================================================
# 
# for i in range(N_v0):
#     v0 = v0_values[i]
#     for j in range(N_t0):
#         theta0 = theta0_values[j]
#         
#         dt = 0.001    
# #        u,T = time_loop_for_cycle(numpy.array([v0, theta0, x0, y0]),dt)
# #        max_distances[i,j]=  u[2]
# #        u,T = time_loop_for_cycle(numpy.array([v0, theta0, x0, y0]),dt)
#         max_distances[i,j]=  distance_theta0_const(v0)
# 
# 
max_dist =  - resbrute[1]    

v0, theta0  = resbrute[0]
vv,tt = resbrute[2]
vv0_values =  vv[:,0]
tt0_values = tt[0,:]


 
print "Maximum distance is %.2f. \
Initial values: v0 = %.2f, t0 = %.2f" % (max_dist,v0,theta0)
 
 
 
dt = 0.001    
u,T = time_loop(numpy.array([v0, theta0, x0, y0]),dt)
 # get the glider's position with respect to the time
x = u[:,2]
y = u[:,3]
v = u[:,0]
 
 # visualization of the path
plt.figure(figsize=(12,6))
plt.grid(True)
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'y', fontsize=18)
plt.title('Glider trajectory, flight time = %.2f , distance = %.2f' % (T,max(x)), fontsize=18)
plt.plot(x,y, 'k-', lw=2);
 
plt.figure(figsize=(12,6))
plt.pcolor(-resbrute[-1])
plt.colorbar()
plt.ylabel(r'$v_0$', fontsize=18)
ylocs,ylabels = plt.yticks()
plt.yticks(ylocs[:-1],tuple('%.2f' % vv0_values[int(i)] for i in ylocs[:-1]))
plt.xlabel(r'$\theta_0$', fontsize=18)
xlocs,xlabels = plt.xticks()
plt.xticks(xlocs[:-1],tuple('%.2f' % degrees(tt0_values[int(i)]) for i in xlocs[:-1]))
plt.suptitle('Maximum distance (in m) for different initial conditions.', fontsize=18)
plt.show()
