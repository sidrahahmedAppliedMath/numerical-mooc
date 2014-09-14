# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 00:47:29 2014

@author: ken
"""
from math import pi
import numpy as np
import matplotlib.pyplot as plt

m_p0 = 100
m_s = 50
g = 9.81
rho = 1.091
r = 0.5
A = pi*r**2
v_e = 325
C_D = 0.15
dt = 0.1


def m_dot_p(t):
    if t<=4.99:
        return 20.
    else:
        return 0
        
def m_p(t):
    if t<=4.99:
        return m_p0 -t*m_dot_p(t)
    else:
        return 0
        
def f(u,t):
    
    v = u[1]
    
    return np.array([v,\
    -g + m_dot_p(t)*v_e/(m_s+m_p(t)) - 0.5*rho*v*abs(v)*A*C_D/(m_s+m_p(t)) \
    ])
    
def euler_step(u, f, dt, t,):
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
    
    return u + dt * f(u,t)
    
u = np.zeros((1, 2))
u[0] = np.array([0.,0.])
T= 0.
n = 0
u = np.append(u,[euler_step(u[n], f, dt,T)],axis = 0)
h = u[n+1,0]
n+=1
T+=dt
h = 1
while h > 0:
    u = np.append(u,[euler_step(u[n], f, dt,T)],axis = 0)
    h = u[n+1,0]
    n+=1
    T+=dt

print u
print m_p(3.2)
time = np.arange(0,T,dt)
y = u[:,0]
vel = u[:,1]
ind_vel_max = np.argmax(vel)
ind_height_max = np.argmax(y)
print "Maximum velocity is %.2f at h = %.2f t = %.2f" % (vel.max(),y[ind_vel_max],time[ind_vel_max])
print "Maximum height is %.2f at t = %.2f with velocity v = %.2f" % (y[ind_height_max],time[ind_height_max],vel[ind_height_max])
print "An impact occured at t = %.2f wit velocity v = %.2f" % (T,vel[-1])
#plt.figure(figsize=(8,6))
#plt.grid(True)
#plt.xlabel(r't', fontsize=18)
#plt.ylabel(r'h', fontsize=18)
##plt.title('Glider trajectory, flight time = %.2f , distance = %.2f' % (T,max(x)), fontsize=18)
#plt.plot(time,y, 'k-', lw=2);