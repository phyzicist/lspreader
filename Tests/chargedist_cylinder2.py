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

def myfunc2(xval, X):
    return np.sum(np.abs(xval - X)) # Sum of distances


    
if __name__ == "__main__":
    print("Hello.")
    
    a = 1.0e-2 # Set a as one cm, 1e4 microns
    rc = 3.0e-6 # Set cylinder radius as 3 microns
    
    # Create the cylinder (rho map)
    l = 5 # Column's line charge density (coulombs/m)
    rho = l/(np.pi*rc**2) # Column's volumetric charge density (coulombs/m^3)
    nx = ny = 100 # Number of points in each dimension
    xgv = np.linspace(-30e-6, 30e-6, nx)
    ygv = np.linspace(-30e-6, 30e-6, ny)
    X, Y = np.meshgrid(xgv, ygv)
    Rc = np.sqrt(X**2 + Y**2) # Radius with respect to cylinder center
    
    Rho = np.zeros_like(Rc)
    Rho[Rc < rc] = rho

    # Convert to linear charge density pixelmap
    dA = np.mean(np.diff(xgv))*np.mean(np.diff(ygv)) # Assuming uniform grid, the area of one cell (m^2)
    L = Rho * dA # Linear charge density (Coulombs/m) for each pixel
    
    # Get an analytical solution, for outside the cylinder
    r0 = 30e-6
    Vabs = (l/(2*np.pi*sc.epsilon_0)) * np.log(r0/Rc)

    # Make a few intermediate-stage plots
    plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    ax.pcolorfast(xgv[[0,-1]]*1e6, ygv[[0,-1]]*1e6, Rho, cmap='viridis')
    plt.title("Rho")
    plt.xlabel("X (um)")
    plt.ylabel("Y (um)")
    plt.axis('equal')
        
    # Numerically solve
    xref = -30.0e-6
    yref = 0.0
    R0 = np.sqrt((X - xref)**2 + (Y - yref)**2) # Radius of each pixel with respect to the designated reference point

    method1 = True
    
    if method1 == True:
        vfunc = np.vectorize(myfunc, excluded=[2,3,4,5])
        V = vfunc(X, Y, X, Y, R0, L)
        
    if method1 == False:
        V = np.zeros_like(Rc)
    
        Xflat = X.flatten()
        Yflat = Y.flatten()
        Vflat = V.flatten()
        Lflat = L.flatten()
        R0flat = R0.flatten()
        Riflat = np.zeros_like(R0flat)
    
        for xval, yval, pot, Xflat2, Yflat2 in np.nditer([X, Y, V, X, Y], op_flags=[['readonly'],['readonly'],['readwrite'],['readonly', 'no_broadcast'], ['readonly','no_broadcast']], flags=['external_loop']): # TODO: Move into external loop for speed?
            #print np.sum(Xflat)
            #print np.sum(Yflat)
            print Xflat2.shape
            Ri = np.sqrt((Xflat2 - xval)**2 + (Yflat2 - yval)**2) # Radius with respect to this point
            print Ri.shape
            ct = (Ri > 0)
            print len(Ri) - ct
                    #print Ri.shape
                    #print xval, yval
                    #print X[0]
                    #print np.sum(Xflat - 1.5e-6)
                    #print(len(Ri) - np.sum(ct))
            #pot[...] = np.sum((Lflat[ct]/(2*np.pi*sc.epsilon_0)) * np.log(R0flat[ct]/Ri[ct])) # Potential at this particular point
                    #pot[...] = np.sqrt(xval**2 + yval**2)        

#    for xval, yval, pot in np.nditer([X, Y, V], op_flags=[['readonly'],['readonly'],['readwrite']], flags=['external_loop']): # TODO: Move into external loop for speed?
#        #print np.sum(Xflat)
#        #print np.sum(Yflat)
#        for i in range(len(Xflat)):
#            Riflat[i] = np.sqrt((Xflat[i] - xval)**2 + (Yflat[i] - yval)**2) # Radius with respect to this point
#        #Ri = np.sqrt((Xflat - xval)**2 + (Yflat - yval)**2) # Radius with respect to this point
#        ct = (Riflat > 0)
#        #print Ri.shape
#        #print xval, yval
#        #print X[0]
#        #print np.sum(Xflat - 1.5e-6)
#        #print(len(Ri) - np.sum(ct))
#        pot[...] = np.sum((Lflat[ct]/(2*np.pi*sc.epsilon_0)) * np.log(R0flat[ct]/Riflat[ct])) # Potential at this particular point
#        #pot[...] = np.sqrt(xval**2 + yval**2)        


#        
#    for j in range(len(xgv)):
#        xval = X[0,j]
#        yval = Y[0,j]
#        print xval, yval
#        Ri = np.sqrt((X - xval)**2 + (Y - yval)**2) # Radius with respect to this point
#        V2[0,j] = np.sum((L/(2*np.pi*sc.epsilon_0)) * np.log(R0/Ri)) # Potential at this particular point

    # Calculate a single one
    xval = -30.0e-6
    yval = -30.0e-6
    Ri = np.sqrt((X - xval)**2 + (Y - yval)**2) # Radius with respect to this point
    ct = Ri > 0
    pot = np.sum((L[ct]/(2*np.pi*sc.epsilon_0)) * np.log( R0[ct] / Ri[ct])) # Potential at this particular point; exclude this point
    print pot/1e11
    
    # Make a few final plots
    plt.figure(3)
    plt.clf()
    plt.plot(X[0]*1e6, Vabs[Rc.shape[0]/2])
    plt.plot(X[0]*1e6, V[Rc.shape[0]/2])
    plt.title("Voltage")
    plt.xlabel("X (um)")
    plt.ylabel("Volts")

    plt.figure(4)
    plt.clf()
    ax = plt.subplot(111)
    ax.pcolorfast(xgv[[0,-1]]*1e6, ygv[[0,-1]]*1e6, V*1e6, cmap='viridis')
    plt.title("Voltage")
    plt.xlabel("X (um)")
    plt.ylabel("Y (um)")
    plt.axis('equal')