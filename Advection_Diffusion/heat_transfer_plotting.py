#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Oct 28 2020

@author: benhills
"""

import numpy as np
import matplotlib.pyplot as plt
from heat_transfer_functions import alpha,surfacePerturbation,harmonicSurface,Robin_T

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

# ---------------------------------------------------------------------------------------------------------------

def numerical_interactive(Ts,zs,ts,maxt,n_profiles):

    plt.figure(figsize=(8,4))

    ax1 = plt.subplot(121)

    # Plot for selected depths
    z_interval = len(zs)//n_profiles
    for i in range(n_profiles):
        idx = i*z_interval
        z = zs[idx]
        T = Ts[1:,idx]
        plt.plot(ts,T,'k',alpha=0.25,label='%.1f'%z)

    plt.ylabel('Temperature ($^\circ$C)');
    plt.xlabel('Years');
    plt.title('Time Series at a Specific Depth');

    ax2 = plt.subplot(122)

    # Plot for selected times
    t_interval = int(len(ts)*(maxt/max(ts)))//n_profiles
    for i in range(n_profiles):
        idx = i*t_interval
        t = ts[idx]
        T = Ts[idx]
        plt.plot(T,zs,'k',alpha=0.25,label='%.2f'%t)

    plt.ylabel('Depth (m)');
    plt.xlabel('Temperature ($^\circ$C)');
    plt.ylim(max(zs),min(zs));
    plt.title('Profile at a Specific Time')

    plt.tight_layout()

    z_idx,t_idx = 0,0
    l1, = ax1.plot(ts,Ts[1:,0],'k')
    p1, = ax1.plot(ts[t_idx],Ts[1+t_idx,0],'k.',mfc='w',ms=10,mew=2)
    l2, = ax2.plot(Ts[0,:],zs,'k')
    p2, = ax2.plot(Ts[0,z_idx],zs[z_idx],'k.',mfc='w',ms=10,mew=2)

    return l1,p1,l2,p2

# ---------------------------------------------------------------------------------------------------------------

def Robin_interactive(Tsurface=-50.,qgeo=0.05,adot=1.,H=1000):

    z,T = Robin_T(Tsurface,qgeo,adot,H)

    plt.figure(figsize=(6,4))
    l1, = plt.plot(T,z,'k')
    plt.xlim(-60,0)
    plt.ylim(0,4000)
    plt.ylabel('Height Above Bed (m)')
    plt.xlabel('Temperature ($^\circ$C)')
    plt.tight_layout()

    return l1
