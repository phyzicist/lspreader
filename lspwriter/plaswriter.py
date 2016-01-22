# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 18:54:17 2016

Make the water column pre-plasma, possibly in a variety of ways in the future.

@author: Scott
"""
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

def simpleColumn(filename, pc_xdims, scale = 1.5, x_crit = -6.07520446, wlen = 800, npoints=100):
    """ Create watercolumn.dat for 1D plasma exponential decay, single scale length.
    Inputs:
        filename: string, path for an output file e.g. 'watercolumn.dat'
        pc_xdims: array-like, (microns) two element array of x limits of plasma, in microns
        scale: number, (microns) the exponential scale length of the pre-plasma
        x_crit: number, (microns) the desired X position of the critical density
        wlen: number, (nanometers) wavelength of light (for calculation of critical density)
        npoints: integer, the number of points to sample for the .dat file
    Outputs:
        Writes a lsp-readable (1D function) file to filename
        X: number, (microns) x values of the pre-plasma
        Y: number, (electrons/cm^3) y values of the pre-plasma
    """
    n_crit = critDens(wlen) # Critical density at this wavelength, in number/cm^3
    n_solid = 1.0e23 # Solid density, in number/cm^3
    
    # Convert dimensions back into LSP units    
    xdims = np.array(pc_xdims)*1e-4 # X dimensions in cm
    scale = scale*1e-4 # Exponential scale length in cm
    x_crit = x_crit * 1e-4 # X position in cm
    xmin = np.min(xdims) 
    xmax = np.max(xdims)

    # Compute the exponential scale length and save the .dat
    X = np.linspace(xmin, xmax, npoints)
    Y = np.zeros(X.shape)
    Y = n_crit * np.exp((X - x_crit)/scale) # Compute a solid density
    Y[Y > n_solid] = n_solid # If density computed exceeds solid density, coerce down to solid density
    Y[0] = Y[-1] = 0 # Set the endpoints to zero
    Xout = -X[::-1] # Funky stuff for save file; reverse axis.
    Yout = Y[::-1]
    np.savetxt(filename,np.array([Xout,Yout]).T,fmt='%.8e',)
    
    X = X*1e4 # Convert back to microns for python output
    return X, Y
    
def columnAnalyze(fn, plot=True):
    """Analyze the 1D watercolumn.dat, by making a plot and printing out critical values at 800 nm
    Inputs:
        fn: string, filename into which to output watercolumn data, e.g. 'watercolumn.dat'
        plot: bool, says whether or not to display a plot of the watercolumn
    Outputs:
        Pops up a plot window (if plot=True) and prints X values of critical, quarter-crit, etc. to the terminal.
    """
    dat = np.loadtxt(fn)
    xgv = -dat[:,0]*1e4 # Convert cm to microns, and reverse
    ne = dat[:,1]
    
    nesold = 1.0 * 10**23 # solid density, electrons per cc
    necrit = 1.742*10**21 # critical density (800 nm), electrons per cc

    if plot:
        plt.plot(xgv, ne)
        plt.axhline(necrit)
        plt.axhline(necrit/4)
        plt.ylim(0,2e21)
        
    _, idx = findNearest(ne, necrit)
    print "Critical density (800nm) is near:", xgv[idx]
    _, idx = findNearest(ne, necrit/4.)
    print "Quarter-critical density is near:", xgv[idx]
    _, idx = findNearest(ne, necrit/(2./3.)**2)
    print "Nine-fourths-critical density is near:", xgv[idx]
    _, idx = findNearest(ne, necrit*4)
    print "Four-critical density is near:", xgv[idx]
    _, idx = findNearest(ne, necrit/10)
    print "Tenth-critical density is near:", xgv[idx]
    
def critDens(wlen):
    """ Calculate the electron critical density (in number/cm^3) for an input vacuum optical wavelength (nm)
    Inputs:
        wlen: Wavelength of light (in nm)
    Outputs:
        n_crit_cm3: Critical electron density (in number/cm^3)
    """
    omega = 2*sc.pi*sc.c/(wlen*1e-9) # Angular frequency of laser in SI units
    n_crit = (sc.epsilon_0*sc.electron_mass/sc.elementary_charge**2)*omega**2 # critical density in SI units (number/m^3)
    n_crit_cm3 = n_crit * 1e-6 # Critical density in number/cm^3
    return n_crit_cm3
    

def findNearest(array, value):
    """
    Get the index for the nearest value to value
    Copied from: http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    Inputs:
        array: 1D array to search for value
        value: number to search for in array
    Outputs:
        nearest value in array, index of that value
    """
    idx = np.abs(array-value).argmin()
    return array[idx], idx