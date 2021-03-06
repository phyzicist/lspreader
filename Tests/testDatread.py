# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 10:33:43 2016

@author: Scott
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b): # Curve fitting exponential function, which must end at 0
    return a * np.exp(-b * x)

def func4(x, a, b, c): # Curve fitting exponential function, without tying to zero
    return a * np.exp(-b * x) + c

def func2(x, a, b): # Curve fitting exponential function
    return a * x + b
 
def func3(x, a, b): # Curve fitting exponential function, which must end at 0
    return a * np.exp(-b * x)

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

cdt = (xgv < -3.5e-4) & (xgv > -10.5e-4)
cline_lim = cline[cdt]
xgv_lim = xgv[cdt]

popt, pcov = curve_fit(func, xgv_lim, cline_lim/np.mean(cline_lim))

## Plot the fitted curve

plt.figure(2)
plt.clf()
#plt.plot(xgv, cline, label="Flash output")
plt.plot(xgv, cline, label="Flash output")
plt.plot(xgv_lim, np.mean(cline_lim)*func(xgv_lim, *popt), 'r--', label="Exponential fit")
plt.ylim(0, 1.5*1.71*1e21/2.65)
plt.title("LSP initial density (Z=0 lineout)")
plt.xlim(-11e-4, -3e-4)
plt.xlabel("LSP X (cm)")
plt.ylabel("LSP triply-ionized H2O electron density (elec/cc)")
plt.hlines([1.71*1e21/2.65, 1.71*1e21/2.65/4], xgv[0], xgv[-1], label='Crit/quarter-crit (fully ionized)', linestyles=['-.'], colors=['k','g'])

plt.legend(loc=2)

## Make the fit to the data

dens_fit = np.zeros(dens.shape)
dens_cline = np.zeros(dens.shape)
for i in range(len(zgv)):
    dens_fit[i,:] = np.mean(cline_lim)*func(xgv, *popt)
    dens_cline[i,:] = cline

# Plot the fitted 2D array
plt.figure(3)
plt.clf()
ax = plt.subplot(111)
im = ax.pcolorfast(X, Z, dens_cline, cmap='viridis', vmax = 1.4e23)
plt.xlabel('X')
plt.ylabel('Z')
plt.colorbar(im)

cline_fit = dens_fit[cix,:] # Lineout through Z = 0

# Replot again
plt.figure(4)
plt.clf()
#plt.plot(xgv, cline, label="Flash output")
plt.plot(xgv, cline, label="FLASH/ LSP central lineout")
plt.plot(xgv, cline_fit, label="Fit/ LSP central lineout")
plt.ylim(0, 1.5*1.71*1e21/2.65)
plt.title("LSP initial density (Z=0 lineout)")
plt.xlim(-11e-4, -3e-4)
plt.xlabel("LSP X (cm)")
plt.ylabel("LSP triply-ionized H2O electron density (elec/cc)")
plt.hlines([1.71*1e21/2.65, 1.71*1e21/2.65/4], xgv[0], xgv[-1], label='Crit/quarter-crit (fully ionized)', linestyles=['-.'], colors=['k','g'])

plt.legend(loc=2)

print "Scale length:", -1/popt[1]*1e4, "microns"

outname = r"C:\Users\Scott\Documents\LSP\Water columns\FLASH sims\f15best_fit.dat"
with open(outname, 'w') as f:
    np.savetxt(f, dens_fit, delimiter = ' ')

outname = r"C:\Users\Scott\Documents\LSP\Water columns\FLASH sims\f15best_cline.dat"
with open(outname, 'w') as f:
    np.savetxt(f, dens_cline, delimiter = ' ')
