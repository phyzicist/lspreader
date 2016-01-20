# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 19:05:18 2016

@author: Scott
"""

import sys
import numpy as np
import lstools as ls
import lspreader2 as rd

# Matplotlib stuff
import matplotlib as mpl
#mpl.use('Agg') # Let's call this at an earlier point, shall we?
import matplotlib.pyplot as plt
from distutils.version import LooseVersion
if LooseVersion(mpl.__version__) < LooseVersion('1.5.0'):    
    # http://stackoverflow.com/questions/11887762/how-to-compare-version-style-strings and 
    print "Matplotlib", mpl.__version__, "might not have colormap 'viridis'. Importing from local colormaps.py."
    import colormaps as cmaps
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    plt.register_cmap(name='inferno', cmap=cmaps.inferno)
    plt.register_cmap(name='magma', cmap=cmaps.magma)
    plt.register_cmap(name='plasma', cmap=cmaps.plasma)

import os
import sftools as sf
import scipy.constants as sc

def addCrit(ax, edens, xgv_um, zgv_um):
    """ Add the critical and quarter-critical density contours over top of current figure axis. Not the best function, but works for now.
    Inputs:    
        ax: axis handle
        edens: electron density NumPy array
        xgv_um: x grid vector, in microns
        zgv_um: z grid vector, in microns
    Outputs:
        Modifies axis 'ax' in place by adding critical density contours to the plot.
    """
    crit = 1.742e21 # critical density (800 nm), electrons per cc
    qcrit = crit/4.0
    
    C = edens
    X, Z = np.meshgrid(xgv_um, zgv_um)
    CS = ax.contour(X, Z, C, [qcrit, crit], linewidths=1.0, colors=('g','w','black'))
    
    fmt = {}
    strs = ['$n_{cr}/4$', '$n_{cr}$', '$n_{s}$']
    for l, s in zip(CS.levels, strs):
        fmt[l] = s
    
    # Label every other level using strings
    ax.clabel(CS, CS.levels, inline_spacing = 10, inline=True, fmt=fmt, fontsize=15, rightside_up=True, manual=[(-20, 13), (-6.2, -10)])
    return CS
    
def addQuiv(ax, vec, xgv, zgv, divsp = 8):
    """ Add a quiver plot of a vector overtop a figure axis.
    Inputs:    
        ax: axis handle
        vec: NumPy array, a vector with three components at every point in 2D (Space x Space x Vector)
        xgv_um: x grid vector, in microns
        zgv_um: z grid vector, in microns
    Outputs:
        Modifies axis 'ax' in place by adding quiver arrows, normalized to themselves, over the plot.
    """
    
    myvecmag = sf.vecMag(vec)
    norm = np.max(myvecmag)
    Vx = vec[::divsp,::divsp,0]/norm
    Vz = vec[::divsp,::divsp,2]/norm
    X, Z = np.meshgrid(xgv[::divsp],zgv[::divsp])
    qu = ax.quiver(X, Z, Vx, Vz, pivot='mid', units='xy', scale_units='xy', scale=0.3, color='white')
    return qu

def mypcolor(C, xgv, zgv, cmin = 0,  cmax = None, title='', tstring = '', clabel = '', fld_id = '', sticker ='', cmap='viridis', edens = np.zeros(0), vec = np.zeros(0)):
    """ A custom-tailored wrapper for pcolorfast(). Somewhat general, meant for any 2D sim colorplot.
    Inputs:
        C: 2D NumPy array, the data to be visualized
        xgv: x dimension grid vector (1D), in microns
        zgv: z dimension grid vector (1D), in microns
        ...
        edens: 2D NumPy array, electron density data for adding critical density contours. If omitted, no contours are added to plot.
    """    
    # xgv, ygv should be in microns
    
    if not cmax:
        cmax = np.max(C)
    fig = plt.figure()
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    xr = [xgv[0],xgv[-1]] # min and max of xgv, needed for pcolorfast
    zr = [zgv[0],zgv[-1]]
    im = ax.pcolorfast(xr, zr, C, cmap=cmap)
    ax.set_xlabel(r'X ($\mu m$)')
    ax.set_ylabel(r'Z ($\mu m$)')
    ax.set_title(title, fontsize=20)
    ax.text(np.min(xgv) + 1, np.max(zgv) - 3, tstring, fontsize=24, color='white')
    ax.text(np.max(xgv) - 6, np.min(zgv) + 2, fld_id, fontsize=44, color='white')
    ax.text(np.min(xgv) + 1, np.min(zgv) + 2, sticker, fontsize=44, color='white')
    cbar = fig.colorbar(im, label=clabel)
    im.set_clim(vmin=cmin, vmax=cmax)
    if len(edens) > 1: # Did the user specify an electron density array?
        addCrit(ax, edens, xgv, zgv) # Highlight the critical density surface
    if len(vec) > 1: # Did the user specify a poynting vector to plot overtop as quivers?
        addQuiv(ax, vec, xgv, zgv)
    return fig

def poyntAnlz(data):
    """ Perform a poynting analysis on a stack of data frames
    Inputs:
        data: the usual data dict; must contain all of the E and B components.
    Outputs:
        Svecmean: array, Full Poynting vector (Space x Space x Vector), averaged across time steps (J/m^2)
        Smagmean: array, Magnitude of Poynting vector (Space x Space), averaged across time steps (J/m^2)
        JEmean: array, Electric field energy density (Space x Space), averaged across time steps (J/m^3)
        JBmean: array, Magnetic field energy density (Space x Space), averaged across time steps (J/m^3)
        Jtotal: number, the total of all electric and magnetic field energy (J/m), averaged across time steps
        """
    # Field energy
    # LSP unit conversion to SI
    tesla = 1e-4 # LSP to SI conversion. 'X Gauss' * tesla = 'Y Tesla'
    vpm = 1e5 # LSP to SI conversion. 'X kV/cm' * vpm = 'Y V/m'
    
    # Poynting vector info: http://hyperphysics.phy-astr.gsu.edu/hbase/waves/emwv.html#c2
    # Energy of electric and magnetic fields: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/engfie.html
    
    Evec = np.stack((data['Ex'], data['Ey'], data['Ez']), axis=-1) * vpm # Electric field vector, in V/m
    Emag = sf.vecMag(Evec)
    Bvec = np.stack((data['Bx'], data['By'], data['Bz']), axis=-1) * tesla # Magnetic field vectors, in Tesla
    Bmag = sf.vecMag(Bvec)
    Svec = (1/sc.mu_0)*np.cross(Evec,Bvec) # Poynting vector, in J/m^2
    Smag = sf.vecMag(Svec)
    
    Svecmean = np.mean(Svec, axis=0) # Mean of the Poynting vector (all components), in J/m^2
    Smagmean = np.mean(Smag, axis=0) # Mean of the Poynting vector's magnitude, in J/m^2
    
    eps = sc.epsilon_0 # Call the dielectric constant equal to vacuum permittivity, which is wrong within a plasma
    mu = sc.mu_0 # Call the magnetic permeability the vacuum permeability, which is wrong within a plasma
    JE = (0.5*eps)*Emag**2 # Electric field energy density, Joule / m^3
    JB = (0.5/mu)*Bmag**2 # Magnetic field energy density, Joule / m^3
    JEmean = np.mean(JE, 0) # Mean Electric field energy density, in J/m^3
    JBmean = np.mean(JB, 0) # Mean Magnetic field energy density, in J/m^3
    
    Asim = (np.max(data['xgv']) - np.min(data['xgv'])) * (np.max(data['zgv']) - np.min(data['zgv'])) * 1e-4 # Area of the sim, in m^2
    Jtotal = np.mean(JEmean + JBmean) * Asim
    print "Simulation total energy: ", Jtotal*1e-3, "mJ / micron."

    return Svecmean, Smagmean, JEmean, JBmean, Jtotal

def plotme(data, outdir='.', alltime=False):
    """ Make Scott's set of custom plots for this batch """
    xgv = data['xgv']*1e4
    zgv = data['zgv']*1e4
    times = data['times']*1e6
    ## CALCULATIONS
    # Mean electron density
    edens = np.mean(data['RhoN10'],0)

    # Mean ion density
    pdens = np.mean(data['RhoN11'],0)
    
    # Mean oxygen ionization state
    old_settings = np.seterr(divide='ignore', invalid='ignore') # Set to ignore divide by zero error. (We will divide by zero where no ions exist)
    ionstate = np.mean((0*data['RhoN1'] + 1*data['RhoN2'] + 2*data['RhoN3'] + 3*data['RhoN4'] + 4*data['RhoN5'] + 5*data['RhoN6'] + 6*data['RhoN7'] + 7*data['RhoN8'])/(data['RhoN1'] + data['RhoN2'] + data['RhoN3'] + data['RhoN4'] + data['RhoN5'] + data['RhoN6'] + data['RhoN7'] + data['RhoN8']), 0)
    ionstate = np.nan_to_num(ionstate)
    np.seterr(**old_settings)  # reset divide by zero error to prior settings
    
    Svec, Smag, JE, JB, Jtot = poyntAnlz(data) # Units are Joules and meters
    JEB_uJum = (JE + JB)*1e-12 # Electromagnetic energy density, uJ/um^3
    Jtot_mJum = Jtot*1e-3 # Total simulation EM energy, in mJ/um

    ## MAKE AND SAVE FIGURES    
    # This timestring computation was copied and pasted from freqanalysis.plotme()
    if alltime: # If this is the all-time rather than time-resolved analysis, put a different label on plot
        maxtime = np.max(times)
        mintime = np.min(times)
        #tstring = r"All time ($\Delta$t =" + "{:.0f}".format(maxtime - mintime) + ' fs)'
        tstring = r"All time (" + "{:.0f}".format(maxtime - mintime) + ' fs total)' # To lay onto the plot
        tlabel = ''.zfill(5) # If this is the "All times" analysis, label the file with "00000.*"
    else: # Otherwise, put the standard "t = XX fs +/- 20 fs" label.
        meantime = np.mean(times)
        tplus = np.max(times) - np.mean(times)
        tstring = 't=' + "{:.1f}".format(meantime) + " fs $\pm$ " + "{:.1f}".format(tplus) + " fs"
        tlabel = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the file label be '00512.*' for t = 51.2342 fs

    pltdir = outdir

    ## Plot 1: Electron density
    C = edens
    sticker = '$e^-$'
    title = 'Electron density'
    clabel = 'Density (number/cc)'
    fld_id = r'$\rho$'
    cmax = 3e21
    fig = mypcolor(C, xgv, zgv, cmin=0,  cmax=cmax, title=title, tstring=tstring, clabel=clabel, fld_id=fld_id, sticker=sticker, edens=edens)
    fig.savefig(os.path.join(sf.subdir(pltdir, 'Electron density'), tlabel + '.png')) # Save into a subdirectory of 'outdir'
    
    ## Plot 2: Proton density
    C = pdens
    sticker = '$p^+$'
    title = 'Proton density'
    clabel = 'Density (number/cc)'
    fld_id = r'$\rho$'
    cmax = 3e21*.67
    fig = mypcolor(C, xgv, zgv, cmin=0,  cmax=cmax, title=title, tstring=tstring, clabel=clabel, fld_id=fld_id, sticker=sticker, edens=edens)
    fig.savefig(os.path.join(sf.subdir(pltdir, 'Proton density'), tlabel + '.png')) # Save into a subdirectory of 'outdir'

    ## Plot 3: Oxygen, mean ionization
    C = ionstate
    sticker = '$O$'
    title = 'Ionization level of oxygen'
    clabel = 'Mean ionization state'
    fld_id = r'+'
    cmin = 0
    cmax = 7
    fig = mypcolor(C, xgv, zgv, cmin=cmin,  cmax=cmax, title=title, tstring=tstring, clabel=clabel, fld_id=fld_id, sticker=sticker, edens=edens, cmap='inferno')
    fig.savefig(os.path.join(sf.subdir(pltdir, 'Oxygen ionization'), tlabel + '.png'))# Save into a subdirectory of 'outdir'
 
    ## Plot 4: Electromagnetic energy density
    C = JEB_uJum
    sticker = r'EM'
    title = 'EM field energy: ' + str(np.round(Jtot_mJum, 2)) + " $mJ/\mu m$"
    clabel = 'Energy density ($nJ/\mu m^3$)'
    fld_id = r'$J$'
    #cmax = 7
    fig = mypcolor(C, xgv, zgv, cmin=0,  cmax=None, title=title, tstring=tstring, clabel=clabel, fld_id=fld_id, sticker=sticker, edens=edens, cmap='inferno', vec=Svec)
    fig.savefig(os.path.join(sf.subdir(pltdir, 'EM energy'), tlabel + '.png'))# Save into a subdirectory of 'outdir'