# Matplotlib stuff (Windows version)
import matplotlib as mpl
import matplotlib.pyplot as plt

from h5stitch2D import *
from freqanalysis import *
import lspreader2 as rd
import numpy as np
import os


def quivme(Svec, Smag, data, pltdir):
    """ Make my custom quiver plot and save to file."""
    times = data['times']*1e6
    xgv = data['xgv']*1e4
    zgv = data['zgv']*1e4

    meantime = np.mean(times)
    tplus = np.max(times) - np.mean(times)
    tstring = 't=' + "{:.1f}".format(meantime) + " fs $\pm$ " + "{:.1f}".format(tplus) + " fs"
    
    ## FIG1: QUIVER PLOT WITH PCOLOR UNDERLAY
    # Quiver documentation: http://matplotlib.org/api/pyplot_api.html
    # Quiver example: http://matplotlib.org/examples/pylab_examples/quiver_demo.html
    fignum = 1
    title = 'Poynting vector, ' + tstring
    vec = Svec
    divsp = 3
    
    figsize = np.array([600,900])/80 # Figure size in inches, assuming 80 dpi
    
    fig = plt.figure(fignum, figsize = figsize)
    plt.clf() # Clear the figure
    ax = plt.subplot(211)
    
    cmin = 0
    cmax = 2124549205619283.2
    
    xr = [xgv[0],xgv[-1]] # min and max of xgv, needed for pcolorfast
    zr = [zgv[0],zgv[-1]]
    C = Smag
    im = ax.pcolorfast(xr, zr, C, cmap='viridis', vmin = cmin, vmax = cmax)
    
    myvecmag = vecMag(vec)
    norm = np.max(myvecmag)
    Vx = vec[::divsp,::divsp,0]/norm
    Vz = vec[::divsp,::divsp,2]/norm
    X, Z = np.meshgrid(xgv[::divsp],zgv[::divsp])
    qu = ax.quiver(X, Z, Vx, Vz, pivot='mid', units='xy', scale_units='xy', scale=1, color='white')
    ax.set_xlabel(r'X ($\mu m$)')
    ax.set_ylabel(r'Z ($\mu m$)')
    ax.set_title(title, fontsize=20)
    #cbar = fig.colorbar(im, label=r'Energy (a.u.)')
    #im.set_clim(vmin=cmin, vmax=cmax)
    
    ## Line plot
    lineout = np.sum(Smag,0)
    ax = plt.subplot(212)
    ax.plot(xgv, lineout)
    ax.set_ylim([0,5e16])
    ax.set_xlabel(r'X ($\mu m$)')
    ax.set_ylabel(r'Poynting magnitude (a.u.)')
    
    fig.set_size_inches(figsize)
    tlabel = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the file label be '00512.*' for t = 51.2342 fs
    fig.savefig(os.path.join(pltdir, tlabel + '.png'))


def vecMag(a, axis=-1):
    """ Gives the magnitude of one (array of) real vector """
    # a: array_like: Components of the vector
    # axis: int, optional: axis that defines the vector. By default, the last axis.
    return np.sqrt(np.sum(np.square(a), axis=axis)) # sqrt(a_x^2 + a_y^2 + ...)
    
p4dir = r'C:\Users\Scott\Documents\temp\2d-nosolid-lowres-'
pltdir = r'C:\Users\Scott\Documents\LSP\Plots\2015-01-13 2d-nosolid-lowres Poynting analysis'

fns = getfnsp4(p4dir)

for stime in range(0,1000,100):
    #stime = 1000
    data = fields2D(fns[stime:stime+100:2])
    
    print data['Ex'].shape
    Evec = np.stack((data['Ex'], data['Ey'], data['Ez']), axis=-1)
    Bvec = np.stack((data['Bx'], data['By'], data['Bz']), axis=-1)
    Svec = np.cross(Evec,Bvec)
    Smag = vecMag(Svec)
    
    Svecmean = np.mean(Svec, axis=0)
    Smagmean = np.mean(Smag, axis=0)
    
    quivme(Svecmean, Smagmean, data, pltdir)
