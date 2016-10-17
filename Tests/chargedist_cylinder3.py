#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
chargedist_cylinder.py: Test problem exploring how quickly I can calculate electric potential of a cylinder of charge

Created by Scott Feister on Mon Oct 10 16:23:43 2016

Will take a 2D profile of a cylinder, but know that it extends in the Z dimension.
Will start with cylinder of large but finite length 2a, and radius R.

Reference for solution: https://www.miniphysics.com/uy1-electric-potential-of-an-infinite-line-charge.html
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

def myfunc(xval, yval, X, Y, R0, L):
    Ri = np.sqrt((X - xval)**2 + (Y - yval)**2) # Radius with respect to this point
    ct = (Ri > 0)
    return np.sum((L[ct]/(2*np.pi*sc.epsilon_0)) * np.log(R0[ct]/Ri[ct])) # Potential at this particular point

def getPot(Xv, Yv, X, Y, Rho, xref=-30.0e-6, yref=0.0):
    """ Get the potential at Xv, Yv points (all inputs and outputs in SI units) """
    # Convert to linear charge density pixelmap
    dA = np.mean(np.diff(X[0,:]))*np.mean(np.diff(Y[:,0])) # Assuming uniform grid, the area of one X/Y/Rho cell (m^2)
    print dA
    L = Rho * dA # Linear charge density (Coulombs/m) for each X/Y/Rho cell

    R0 = np.sqrt((X - xref)**2 + (Y - yref)**2) # Radius of each pixel with respect to the designated reference point

    # Numerically solve
    vfunc = np.vectorize(myfunc, excluded=[2,3,4,5])
    Vv = vfunc(Xv, Yv, X, Y, R0, L)

    return Vv

    
if __name__ == "__main__":
    print("Hello.")
    
    # Create the cylinder (rho map)
    rc = 3.0e-6 # Set cylinder radius as 3 microns
    l = 5 # Column's line charge density (coulombs/m)
    rho = l/(np.pi*rc**2) # Column's volumetric charge density (coulombs/m^3)
    nx = ny = 800 # Number of points in each dimension
    xgv = np.linspace(-30e-6, 30e-6, nx)
    ygv = np.linspace(-30e-6, 30e-6, ny)
    X, Y = np.meshgrid(xgv, ygv)
    Rc = np.sqrt(X**2 + Y**2) # Radius with respect to cylinder center
    
    Rho = np.zeros_like(Rc)
    Rho[Rc < rc] = rho

    # Solve for the potential
    Xv, Yv = np.meshgrid(xgv, [0.0])
    
    Vv = getPot(Xv, Yv, X, Y, Rho)

    # Get an analytical solution, for outside the cylinder
    r0 = 30.0e-6
    Vabs = (l/(2*np.pi*sc.epsilon_0)) * np.log(r0/Rc)

    # Make a few plots
    plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    ax.pcolorfast(xgv[[0,-1]]*1e6, ygv[[0,-1]]*1e6, Rho, cmap='viridis')
    plt.title("Rho")
    plt.xlabel("X (um)")
    plt.ylabel("Y (um)")
    plt.axis('equal')
    
    plt.figure(2)
    plt.clf()
    plt.plot(X[0]*1e6, Vabs[Vabs.shape[0]/2])
    plt.plot(Xv[0]*1e6, Vv[Vv.shape[0]/2])
    plt.title("Voltage")
    plt.xlabel("X (um)")
    plt.ylabel("Volts")

    plt.figure(3)
    plt.clf()
    ax = plt.subplot(111)
    ax.pcolorfast(xgv[[0,-1]]*1e6, ygv[[0,-1]]*1e6, Vv, cmap='viridis')
    plt.title("Voltage")
    plt.xlabel("X (um)")
    plt.ylabel("Y (um)")
    plt.axis('equal')