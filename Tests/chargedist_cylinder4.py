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
import sftools as sf
import time

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

def getPot2(Xv, Yv, X, Y, Rho, xref=-30.0e-6, yref=0.0):
    """ Get the potential at Xv, Yv points (all inputs and outputs in SI units) 
    BEST ALGO so far.
    """
    # Convert to linear charge density pixelmap
    dA = np.mean(np.diff(X[0,:]))*np.mean(np.diff(Y[:,0])) # Assuming uniform grid, the area of one X/Y/Rho cell (m^2)
    L = Rho * dA # Linear charge density (Coulombs/m) for each X/Y/Rho cell

    vshape = Xv.shape
    R0 = np.sqrt((X - xref)**2 + (Y - yref)**2) # Radius of each pixel with respect to the designated reference point

    # Static arrays
    X = X.flatten()
    Y = Y.flatten()
    L = L.flatten()
    R0 = R0.flatten()
    
    # Temporary array
    Ri = np.zeros_like(R0)

    # Arrays to be chunked over
    Xv = Xv.flatten()
    Yv = Yv.flatten()
    Vv = np.zeros_like(Xv)

    # Now do the magic
    for i in range(len(Xv)):
        Ri = np.sqrt((X - Xv[i])**2 + (Y - Yv[i])**2) # Radius with respect to this point
        ct = (Ri > 0)
        Vv[i] = np.sum((L[ct]/(2*np.pi*sc.epsilon_0)) * np.log(R0[ct]/Ri[ct])) # Potential at this particular point

    Vv = np.reshape(Vv, vshape)
    return Vv

def getPot4(Xv, Yv, X, Y, Rho, xref=-30.0e-6, yref=0.0, cksz=5):
    """ Get the potential at Xv, Yv points (all inputs and outputs in SI units). Try chunking. 
    
    cksz = 200 # Size of the chunks for memory
    Insanely memory intensive algorithm, and not faster. Ditch. BAD ALGO.
    """
    # Convert to linear charge density pixelmap
    dA = np.mean(np.diff(X[0,:]))*np.mean(np.diff(Y[:,0])) # Assuming uniform grid, the area of one X/Y/Rho cell (m^2)
    L = Rho * dA # Linear charge density (Coulombs/m) for each X/Y/Rho cell

    vshape = Xv.shape
    R0 = np.sqrt((X - xref)**2 + (Y - yref)**2) # Radius of each pixel with respect to the designated reference point

    # Static arrays
    X = np.tile(X.flatten(), (cksz, 1)).T
    Y = np.tile(Y.flatten(), (cksz, 1)).T
    L = np.tile(L.flatten(), (cksz, 1)).T
    R0 = np.tile(R0.flatten(), (cksz, 1)).T
    
    print("R0 shape: " + str(R0.shape))
    # Temporary array
    Ri = np.zeros_like(R0)

    # Arrays to be chunked over
    Xv = Xv.flatten()
    Yv = Yv.flatten()
    Vv = np.zeros_like(Xv)

    # Now do the magic
    ix_all = range(len(Xv))
    ix_chunks = sf.chunkFx(ix_all, cksz)[:-1] # Kill off the odd-sized chunk, for now
    
    for ixs in ix_chunks:
        Ri = np.sqrt((X - Xv[ixs])**2 + (Y - Yv[ixs])**2) # Radius with respect to this point
        ct = (Ri > 0)
        Vv[ixs] = np.sum((L[ct]/(2*np.pi*sc.epsilon_0)) * np.log(R0[ct]/Ri[ct])) # Potential at this particular point

    Vv = np.reshape(Vv, vshape)
    return Vv
    
if __name__ == "__main__":
    print("Hello.")
    
    # Create the cylinder (rho map)
    rc = 6.0e-6 # Set cylinder radius as 3 microns
    l = 5 # Column's line charge density (coulombs/m)
    rho = l/(np.pi*rc**2) # Column's volumetric charge density (coulombs/m^3)
    nx = ny = 100 # Number of points in each dimension
    xgv = np.linspace(-30e-6, 30e-6, nx)
    ygv = np.linspace(-30e-6, 30e-6, ny)
    X, Y = np.meshgrid(xgv, ygv)
    Rc = np.sqrt(X**2 + Y**2) # Radius with respect to cylinder center
    
    Rho = np.zeros_like(Rc)
    Rho[Rc < rc] = rho

    # Solve for the potential
    #Xv, Yv = np.meshgrid(xgv, [0.0])
    Xv, Yv = X, Y
    
#    t0 = time.time()
#    Vv = getPot4(Xv, Yv, X, Y, Rho)
#    t1 = time.time()
#    print("getPot4: " + str(t1 - t0) + " secs")

    t0 = time.time()
    Vv = getPot2(Xv, Yv, X, Y, Rho)
    t1 = time.time()
    print("getPot2: " + str(t1 - t0) + " secs")
    
    # Gradient calculation
    dxv = np.mean(np.diff(Xv[0,:]))    
    dyv = np.mean(np.diff(Yv[:,0]))
    [Eyv, Exv] = np.gradient(-Vv, dyv, dxv) # Is this order swapped???
    
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

    plt.figure(4)
    plt.clf()
    skip = 5
    Q = plt.quiver(Exv[::skip,::skip], Eyv[::skip,::skip])
    plt.title("Electric field")
    plt.xlabel("X (um)")
    plt.ylabel("Y (um)")
    plt.axis('equal')
