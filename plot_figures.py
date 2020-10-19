#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Oct 18 2020

@author: benhills
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_perturbation_figure(T0,dT,t,zs,alpha,surfacePerturbation):
    # Run the function
    T = surfacePerturbation(T0,dT,t,zs,alpha)

    # Plot
    plt.figure(figsize=(4,4))
    plt.plot(T,zs,'k')
    plt.ylim(max(zs),min(zs))
    plt.ylabel('Depth (m)');
    plt.xlabel('Ice Temperature ($^\circ$C)');


def plot_harmonic_figures(zs,ts,T0,Ta,omega,n_profiles,harmonicSurface):

    plt.figure(figsize=(4,3))
    plt.title('Surface Temperature')
    plt.plot(ts,T0+Ta*np.sin(2.*np.pi*omega*ts),'k');
    plt.ylabel('$^\circ$C')
    plt.xlabel('Years');
    plt.tight_layout()

    plt.figure(figsize=(8,4))
    plt.subplot(121)
    # Run and plot for all the times
    for t in ts[::len(ts)//n_profiles]:
        # Run the function
        T = harmonicSurface(T0,Ta,t,zs,omega)
        # Plot
        plt.plot(T,zs,'k',alpha=0.25+0.75*(t/max(ts)),label='%.2f'%(t))

    plt.legend(title='Years',bbox_to_anchor=(1, 1))
    plt.ylabel('Depth (m)');
    plt.xlabel('Ice Temperature ($^\circ$C)');
    plt.ylim(max(zs),min(zs));
    plt.title('Profile at a Specific Time')

    plt.subplot(122)
    # Plot the time series for all depths
    for z in zs[::len(zs)//n_profiles]:
        T = harmonicSurface(T0,Ta,ts,z,omega)
        # Plot
        plt.plot(ts,T,'k',alpha=0.25+0.75*float(z)/max(zs),label='%.0f'%z)

    plt.legend(title='Depths (m)',bbox_to_anchor=(1, 1))
    plt.ylabel('Ice Temperature ($^\circ$C)');
    plt.xlabel('Years');
    plt.title('Time Series at a Specific Depth');
    plt.tight_layout()
