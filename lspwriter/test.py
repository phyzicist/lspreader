# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 14:33:23 2016

@author: Scott
"""

from plaswriter import radialPlume, parabCup
import matplotlib.pyplot as plt
import numpy as np

scale = 2.8#um # Exponential scale length of pre-plasma
pc_xdims = (-30, 0)#um # X limits to simulation particle creation space, in microns
pc_zdims = (-15, 15)#um # Z limits to simulation particle creation space, in microns

plumedat = r'C:\Users\Scott\Documents\LSP\LSPwrite submissions\TEST\2Dplume.dat'

#D, xgv, zgv = radialPlume(plumedat, pc_xdims, pc_zdims, scale = scale, nxpoints = 100, nzpoints = 400)
D, xgv, zgv = parabCup(plumedat, pc_xdims, pc_zdims, parabfoc = 1.5, scale = scale, nxpoints = 100, nzpoints = 400)

X, Z = np.meshgrid(xgv, zgv)

plt.figure(1)
plt.clf()
plt.contourf(X, Z, np.log10(D), [np.log10(1.74e21/4), np.log10(1.74e21), 22.9, 24], cmap='viridis')
plt.colorbar()
plt.axis('equal')

plt.figure(2)
plt.clf()
plt.contourf(X, Z, np.log10(D), cmap='viridis')
plt.colorbar()
plt.axis('equal')