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

### Create model matrix for a semi-implicit Crank-Nicholson scheme ###

from scipy import sparse
def sparseMatrix(dt,zs,u,qgeo,conductivity=2.1):
    # calculate the diffusion and advection terms
    dz = np.mean(np.gradient(zs))
    diff = alpha*(dt/(dz**2.))/2.
    adv = u*dt/(4.*dz)

    # Write the sparse matrices for left and rhs of equation, A*Tnew=B*Tlast+S
    N = len(zs)
    A = sparse.lil_matrix((N, N))           # Create a sparse Matrix
    A.setdiag(1.+2.*diff*np.ones(N))        # Set the diagonals
    A.setdiag(-(diff-adv)*np.ones(N),k=1)
    A.setdiag(-(diff+adv)*np.ones(N),k=-1)
    B = sparse.lil_matrix((N, N))
    B.setdiag(1.-2.*diff*np.ones(N))
    B.setdiag((diff-adv)*np.ones(N),k=1)
    B.setdiag((diff+adv)*np.ones(N),k=-1)

    # Boundary Conditions
    A[-1,-1] = 1+2.*diff    # Neumann at bed
    A[-1,-2] = -2.*diff
    B[-1,-1] = 1-2.*diff
    B[-1,-2] = 2.*diff
    A[0,:] = 0.             # Dirichlet at surface
    A[0,0] = 1.
    B[0,:] = 0.
    B[0,0] = 1.

    # geothermal source
    S = np.zeros(N)
    S[-1] = -dt*qgeo/conductivity*(2./dz)*alpha

    return A.tocsr(),B.tocsr(),S

### Numerical Model ###

from scipy.sparse.linalg import spsolve
def numericalModel(zs,ts,dt,u=0,BC_upper=[0.],BC_lower=0.,IC=[None]):
    # Initial condition
    if np.all(IC) == None:
        T = BC_upper[0]*np.ones_like(zs)
    else:
        T = IC
    if len(BC_upper) == 1:
        BC_upper = BC_upper[0]*np.ones_like(ts)
    # Set up the matrices
    A,B,S = sparseMatrix(dt*spy,zs,u/spy,BC_lower)
    # loop through times
    T_out = np.array([T])
    for i in range(len(ts)):
        T[0] = BC_upper[i]
        rhs = B*T+S
        T = spsolve(A,rhs)
        T_out = np.append(T_out,[T],axis=0)
    return T_out

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
