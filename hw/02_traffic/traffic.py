# -*- coding: utf-8 -*-
"""
Created on Sat Oct  4 23:56:49 2014

@author: ken
"""

import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

v_max = 80. #Part A
#v_max = 136.
L= 11.
rho_max = 250.
nx =51
dt = .001


dx =L/(nx-1)

#initial conditions

x = numpy.linspace(0,L,nx)
rho0 = numpy.ones(nx)*10 # Part A
#rho0 = numpy.ones(nx)*20 # Part B
rho0[10:20] = 50

def velocity(rho):
    return v_max*(1 - rho/rho_max)
    
def Jacobian(rho):
    return v_max*(1 - 2*rho/rho_max)


def Flux(rho):
    return rho*velocity(rho)
    
t = 6./60    
nt = int(t/dt+1)
rho = rho0.copy()    
for n in range(1, nt):  
    rhon = rho.copy() 
    rho[1:] = rhon[1:]-0.5*(Jacobian(rhon)[1:]+Jacobian(rhon)[0:-1])*dt/dx*(rhon[1:]-rhon[0:-1]) 
    rho[0] = rho0[0]

#rho = rho0.copy()    
#for n in range(1, nt):  
#    rhon = rho.copy() 
#    Fluxn = Flux(rhon)
#    rho[1:] = rhon[1:]-dt/dx*(Fluxn[1:]-Fluxn[0:-1]) 
#    rho[0] = rho0[0]
    
vel =  velocity(rho)
ave_vel = numpy.average(vel)
print "Minimum initial velocity : {:.02f}".format(min(velocity(rho0))/3.6) 
print "Minimum velocity : {:.02f}".format(min(vel)/3.6) 
print "Average velocity : {:.02f}".format(ave_vel/3.6) 
plt.plot(x,velocity(rho))