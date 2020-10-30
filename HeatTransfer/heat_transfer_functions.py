#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Oct 18 2020

@author: benhills
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import erf
from constants import Constants

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

def Robin_T(Ts,qgeo,H,adot,nz=101,
        const=Constants(),melt=True,verbose=True):

    """
    Analytic ice temperature model from Robin (1955)

    Assumptions:
        1) no horizontal advection
        2) vertical advection is linear with depth
        3) firn column is treated as equivalent thickness of ice
        4) If base is warmer than the melting temperature recalculate with new basal gradient
        5) no strain heating

    Parameters
    ----------
    Ts:     float,  Surface Temperature (C)
    qgeo:   float,  Geothermal flux (W/m2)
    H:      float,  Ice thickness (m)
    adot:   float,  Accumulation rate (m/yr)
    nz:     int,    Number of layers in the ice column
    const:  class,  Constants
    melt:   bool,   Choice to allow melting, when true the bed temperature
                    is locked at the pressure melting point and melt rates
                    are calculated

    Output
    ----------
    z:      1-D array,  Discretized height above bed through the ice column
    T:      1-D array,  Analytic solution for ice temperature
    """

    z = np.linspace(0,H,nz)
    adot/=const.spy
    q2 = adot/(2*(const.k/(const.rho*const.Cp))*H)
    Tb_grad = -qgeo/const.k
    f = lambda z : np.exp(-(z**2.)*q2)
    TTb = Tb_grad*np.array([quad(f,0,zi)[0] for zi in z])
    dTs = Ts - TTb[-1]
    T = TTb + dTs
    # recalculate if basal temperature is above melting (see van der Veen pg 148)
    Tm = const.beta*const.rho*const.g*H
    if melt and T[0] > Tm:
        Tb_grad = -2.*np.sqrt(q2)*(Tm-Ts)/np.sqrt(np.pi)*(np.sqrt(erf(adot*H*const.rho*const.Cp/(2.*const.k)))**(-1))
        TTb = Tb_grad*np.array([quad(f,0,zi)[0] for zi in z])
        dTs = Ts - TTb[-1]
        T = TTb + dTs
    return z,T
