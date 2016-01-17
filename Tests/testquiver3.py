# Matplotlib stuff (Windows version)
import matplotlib as mpl
import matplotlib.pyplot as plt

from h5stitch2D import *
from freqanalysis import *
import lspreader2 as rd
import numpy as np
import os
import scipy.constants as sc # Scientific constants e.g. speed of light, vacuum permittivity (in SI units)

def quivme(Svec, Smag, data, pltdir, fignum):
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
    title = 'Poynting vector, ' + tstring
    vec = Svec
    divsp = 3
    
    figsize = np.array([600,900])/80 # Figure size in inches, assuming 80 dpi
    
    fig = plt.figure(fignum, figsize = figsize)
    plt.clf() # Clear the figure
    ax = plt.subplot(211)
    

    
    xr = [xgv[0],xgv[-1]] # min and max of xgv, needed for pcolorfast
    zr = [zgv[0],zgv[-1]]
    C = Smag
    cmin = 0
    cmax = np.max(C)
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
    #ax.set_ylim([0,5e16])
    ax.set_xlabel(r'X ($\mu m$)')
    ax.set_ylabel(r'Poynting magnitude (a.u.)')
    
    fig.set_size_inches(figsize)
    tlabel = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the file label be '00512.*' for t = 51.2342 fs
    #fig.savefig(os.path.join(pltdir, tlabel + '.png')) ## TURNED OFF PLOT MAKING!

def vecMag(a, axis=-1):
    """ Gives the magnitude of one (array of) real or complex vector """
    # a: array_like: Components of the vector
    # axis: int, optional: axis that defines the vector. By default, the last axis.
    return np.sqrt(np.sum(np.multiply(a,np.conj(a)), axis=axis)) # sqrt((a_x a_x*)+ (a_y a_y*) + ...)
    
p4dir = r'C:\Users\Scott\Documents\temp\2d-nosolid-lowres-'
pltdir = r'C:\Users\Scott\Documents\LSP\Plots\2015-01-13 2d-nosolid-lowres Poynting analysis'

#fns = getfnsp4(p4dir)
fns = [r'C:\Users\Scott\Documents\temp\Jan 14 tests\5um_Foc-14 comparisons\Lres\flds250.p4']
#fns = [r'C:\Users\Scott\Documents\temp\Jan 14 tests\5um_Foc-14 comparisons\Hres\flds625.p4']
data = fields2D(fns[::1])

# LSP unit conversion to SI
tesla = 1e-4 # LSP to SI conversion. 'X Gauss' * tesla = 'Y Tesla'
vpm = 1e5 # LSP to SI conversion. 'X kV/cm' * vpm = 'Y V/m'

# Poynting vector info: http://hyperphysics.phy-astr.gsu.edu/hbase/waves/emwv.html#c2
# Energy of electric and magnetic fields: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/engfie.html

print data['Ex'].shape
Evec = np.stack((data['Ex'], data['Ey'], data['Ez']), axis=-1) * vpm # Electric field vector, in V/m
Emag = vecMag(Evec)
Bvec = np.stack((data['Bx'], data['By'], data['Bz']), axis=-1) * tesla # Magnetic field vectors, in Tesla
Bmag = vecMag(Bvec)
Svec = (1/sc.mu_0)*np.cross(Evec,Bvec) # Poynting vector, in J/m^2
Smag = vecMag(Svec)


Svecmean = np.mean(Svec, axis=0)
Smagmean = np.mean(Smag, axis=0)

eps = sc.epsilon_0 # Call the dielectric constant equal to vacuum permittivity, which is wrong within a plasma
mu = sc.mu_0 # Call the magnetic permeability the vacuum permeability, which is wrong within a plasma
JEmag = (0.5*eps)*Emag**2 # Electric field energy density, Joule / m^3
JBmag = (0.5/mu)*Bmag**2 # Magnetic field energy density, Joule / m^3

Asim = (np.max(data['xgv']) - np.min(data['xgv'])) * (np.max(data['zgv']) - np.min(data['zgv'])) * 1e-4 # Area of the sim, in m^2
Vsim = Asim # Have to stuff in the height of the sim somewhere
Jtotal = np.mean(JBmag[0,:,:] + JEmag[0,:,:]) * Asim # Joules / m (since we can't get rid of that last dimension)
print "Simulation total energy: ", Jtotal*1e-3, "mJ / micron."
quivme(Svecmean, Smagmean, data, pltdir, 1)

title = "Energy density"
fignum = 2
C = JBmag[0,:,:] + JEmag[0,:,:]
fig = plt.figure(fignum)
plt.clf() # Clear the figure
cmin = np.min(C)
cmax = np.max(C)
ax = plt.subplot(111)
xr = [np.min(xgv),np.max(xgv)] # min and max of xgv, needed for pcolorfast
zr = [np.min(zgv),np.max(zgv)]
im = ax.pcolorfast(xr,zr,C, cmap='viridis')
ax.set_xlabel(r'X ($\mu m$)')
ax.set_ylabel(r'Z ($\mu m$)')
ax.set_title(title, fontsize=20)
cbar = fig.colorbar(im, label=r'Energy density (J/m$^3$)')
im.set_clim(vmin=cmin, vmax=cmax)
    
