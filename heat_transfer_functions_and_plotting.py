#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Oct 18 2020

@author: benhills
"""

import numpy as np
import matplotlib.pyplot as plt

# Dimensional Constant
spy  = 60.*60.*24.*365.24    # seconds per year (s yr-1)
alpha = 1.09e-6                  # Thermal Diffusivity of Ice (m2 s-1), if you change this the actual value is 1.09e-6

# ---------------------------------------------------------------------------------------------------------------

# Function for the surface perturbation
def surfacePerturbation(T0,dT,t,z,alpha):
    t *= spy
    if alpha > .1:
        alpha /= spy
    from scipy.special import erfc
    # equation 2.4-7 from Carslaw and Jaegar (1959)
    T = T0 + dT*erfc(abs(z)/(2.*np.sqrt(alpha*t)))
    return T

# Surface Temperature Harmonic Function
def harmonicSurface(T0,Ta,t,z,omega,alpha=alpha):
    alpha *= spy
    # equation 2.6-8 from Carslaw and Jaegar (1959)
    T = T0 + Ta * np.exp(-z*np.sqrt(np.pi*omega/(alpha))) * np.sin((2.*np.pi*omega*t)-z*np.sqrt(np.pi*omega/(alpha)))
    return T

# ---------------------------------------------------------------------------------------------------------------

def perturbation_interactive(T0=0.,dT=1.,t=1.,zs=np.linspace(0,1000,100),
                        alpha=alpha,surfacePerturbation=surfacePerturbation):

    # Run the function
    T = surfacePerturbation(T0,dT,t,zs,alpha)

    # Plot
    plt.figure(figsize=(4,4))
    l1, = plt.plot(T,zs,'k')
    plt.ylim(max(zs),min(zs))
    plt.xlim(-1,10)
    plt.ylabel('Depth (m)');
    plt.xlabel('Ice Temperature ($^\circ$C)');

    plt.tight_layout()

    return l1

# ---------------------------------------------------------------------------------------------------------------

def harmonic_interactive(zs=np.linspace(0,20,100),ts=np.linspace(0,1,100),T0=0.,
                    Ta=1.,omega=1.,n_profiles=10,harmonicSurface=harmonicSurface):

    plt.figure(figsize=(8,4))

    ax1 = plt.subplot(121)
    # Plot the time series for all depths
    for z in zs[::len(zs)//n_profiles]:
        T = harmonicSurface(T0,Ta,ts,z,omega)
        # Plot
        plt.plot(ts,T,'k',alpha=0.1,label='%.0f'%z)

    plt.ylabel('Ice Temperature ($^\circ$C)');
    plt.xlabel('Years');
    plt.xlim(0,1/omega)
    plt.title('Time Series at a Specific Depth');

    ax2 = plt.subplot(122)
    # Run and plot for all the times
    for t in ts[::len(ts)//n_profiles]:
        # Run the function
        T = harmonicSurface(T0,Ta,t,zs,omega)
        # Plot
        plt.plot(T,zs,'k',alpha=0.1,label='%.2f'%(t))

    plt.ylabel('Depth (m)');
    plt.xlabel('Ice Temperature ($^\circ$C)');
    plt.ylim(max(zs),min(zs));
    plt.xlim(-1.1,1.1)
    plt.title('Profile at a Specific Time')

    plt.tight_layout()

    z_init,t_init = 0,0
    T = harmonicSurface(T0,Ta,ts,z_init,omega)
    l1, = ax1.plot(ts,T,'k')
    p1, = ax1.plot(t_init,T[np.argmin(abs(ts-t_init))],'k.',mfc='w',ms=10,mew=2)
    T = harmonicSurface(T0,Ta,t_init,zs,omega)
    l2, = ax2.plot(T,zs,'k')
    p2, = ax2.plot(T[np.argmin(abs(zs-z_init))],z_init,'k.',mfc='w',ms=10,mew=2)

    return l1,p1,l2,p2
