# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 10:33:43 2016

@author: Scott
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import signal

def func(x, a, b): # Curve fitting exponential function, which must end at 0
    return a * np.exp(-b * x)

#def func(x, a, b, c, d): # Curve fitting third-order polynomial function
#    return a + b * x + c * x**2 + d * x**3

fn = r"C:\Users\Scott\Documents\LSP\Water columns\FLASH sims\f15best.dat"

xgv = np.genfromtxt(fn, skip_header=4, skip_footer=404)
zgv = np.genfromtxt(fn, skip_header=5, skip_footer=403)
dens = np.genfromtxt(fn, skip_header=6)
print xgv.shape
print zgv.shape
print dens.shape

X, Z = np.meshgrid(np.concatenate((xgv, xgv[-1:])), np.concatenate((zgv, zgv[-1:])))

## Plot the initial setup

plt.figure(1)
plt.clf()
ax = plt.subplot(111)
im = ax.pcolorfast(X, Z, dens, cmap='viridis', vmax = 1.4e23)
plt.xlabel('X')
plt.ylabel('Z')
plt.colorbar(im)

# Extract a lineout
cix = np.int(np.round(len(zgv)/2)) # Index of Z = 0
cline = dens[cix,:] # Lineout through Z = 0


## Make the fit to the data
cdt = (xgv < -8e-4) & (xgv > -30e-4)
cline_lim = cline[cdt]
popt, pcov = curve_fit(func, xgv[cdt], cline_lim/np.mean(cline_lim))
ep1line = np.mean(cline_lim)*func(xgv, *popt)

## Make the fit to the data
cdt = (xgv < -8e-4) & (xgv > -30e-4)
cline_lim = cline[cdt]
popt, pcov = curve_fit(func, xgv[cdt], cline_lim/np.mean(cline_lim))
ep1line = np.mean(cline_lim)*func(xgv, *popt)
sc1 = -1/popt[1]*1e4
print "Fit 1 scale length:", sc1, "microns"

## Make the fit to the data
cdt = (xgv > -5e-4) & (xgv < -1e-4)
cline_lim = cline[cdt]
popt, pcov = curve_fit(func, xgv[cdt], cline_lim/np.mean(cline_lim))
ep2line = np.mean(cline_lim)*func(xgv, *popt)
sc2 = -1/popt[1]*1e4
print "Fit 2 scale length:", sc2, "microns"

jx = np.argmax(xgv > -10e-4)
jx2 = np.argmax(xgv > -26e-4)
dens_fit = np.zeros(dens.shape)
for i in range(len(zgv)):
    dens_fit[i,:] = cline
    dens_fit[i,:jx] = ep2line[:jx]
    slope = (dens_fit[i,jx2] - 0)/(xgv[jx2]  - xgv[0])
    dens_fit[i,:jx2] = slope * (xgv[:jx2] - xgv[0])
    dens_fit[i,:] = signal.savgol_filter(dens_fit[i,:], 25, 3)


# Replot again
plt.figure(4)
plt.clf()
#plt.plot(xgv, cline, label="Flash output")
plt.plot(xgv*1e4, cline, label="FLASH density, central linout")
plt.plot(xgv*1e4, dens_fit[cix,:], 'g', label="Modified with " + str(np.round(sc2,1)) + " $\mu$m scale length decay")
#plt.plot(xgv*1e4, ep1line, 'r--', label=str(np.round(sc1,1)) + " $\mu$m scale length fit")

plt.ylim(0, 1*1.71*1e21/2.65)
plt.title("LSP initial density (Z=0 lineout)")
plt.xlim(-30, 0)
plt.xlabel("LSP X ($\mu$m)")
plt.ylabel("LSP triply-ionized H2O electron density (elec/cc)")
plt.hlines([1.71*1e21/2.65, 1.71*1e21/2.65/4], xgv[0]*1e4, xgv[-1]*1e4, label='Crit/quarter-crit (fully ionized)', linestyles=['-.'], colors=['k','g'])

plt.legend(loc=2)


outname = r"C:\Users\Scott\Documents\LSP\Water columns\FLASH sims\f15vacdecay.dat"
with open(outname, 'w') as f:
    np.savetxt(f, dens_fit, delimiter = ' ')