# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 15:20:21 2016

@author: Scott
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def exp1(x, a, l, x0): # Exponential function
    # Rising exponential, scale length l, with y(x0) = a
    return a * np.exp((x - x0)/l)


def twoscale(x, l1 = 1.0, l2 = 20.0, n2 = 1.6*0.285e21/2.65, xc = -4.6):
    """
    x: 1D array of X values on which to calculate
    xc: 1D X value of critical density
    n2: Density at which to switch scales
    l1: Steep exponential scale length
    l2: Slow exponential scale length
    
    TODO:
        * Do I need to round off the edges? Unclear
    """
    ns = 1.0e23/2.65 # Solid density (accounts for scaling in LSP load)
    nc = 1.71e21/2.65 # Critical density (accounts for scaling in LSP load)
    
    # Calculate the new scale
    y = exp1(x, nc, l1, xc)
    x2 = x[np.argmax(y > n2)] # X value to switch to slow scale length
    y2 = exp1(x, n2, l2, x2)
    y[y>ns] = ns
    y[x < x2] = y2[x < x2]
    return y

def threescale(x, l0 = 0.2, l1 = 1.0, l2 = 20.0, n0 = 9e21/2.65, n2 = 1.6*0.285e21/2.65, xc = -4.6):
    """
    x: 1D array of X values on which to calculate
    xc: 1D X value of critical density
    n2: Density at which to switch scales
    l1: Steep exponential scale length
    l2: Slow exponential scale length
    
    TODO:
        * Do I need to round off the edges? Unclear
    """
    ns = 3*1.0e23/2.65 # Solid density / shock density (accounts for scaling in LSP load)
    nc = 1.71e21/2.65 # Critical density (accounts for scaling in LSP load)
    
    # Calculate the new scale
    y = exp1(x, nc, l1, xc)
    
    x2 = x[np.argmax(y > n2)] # X value to switch to slow scale length
    y2 = exp1(x, n2, l2, x2)

    x0 = x[np.argmax(y > n0)] # X value to switch to steep scale length
    y0 = exp1(x, n0, l0, x0)


    y[x < x2] = y2[x < x2]
    y[x > x0] = y0[x > x0]
    y[y > ns] = ns

    return y


def write1DLSP(fname, x, y):
    # x should be in microns units
    # Convert dimensions back into LSP units
    x = x * 1e-4
    
    with open(fname, 'w') as f:
        f.write('# 1D LSP dat file (function type 30)\n')
        f.write('# X (cm), Y\n')
        np.savetxt(f,np.array([-x[::-1],y[::-1]]).T,fmt='%.8e',)

def write2DLSP(fname, xgv, zgv, D):
    # xgv and zgv should be in LSP units
    # D is the 2D array to write
    # Convert dimensions back into LSP units
    xgv = xgv * 1e-4
    zgv = zgv * 1e-4

    with open(fname, 'w') as f:
        f.write('# 2D LSP dat file (function type 40)\n')
        f.write('# Dimensions X, Z, Density. See page 147 of 200 in LSP manual, or search "type 40", for specification format.\n')
        f.write(str(int(xgv.shape[0])) + ' ' + str(int(zgv.shape[0])) + '\n')
        np.savetxt(f, xgv[None], delimiter = ' ') # xgv[None] changes the dimensions of xgv from (100,) to (1, 100); this lets us save it as a row rather than column
        np.savetxt(f, zgv[None], delimiter = ' ')
        np.savetxt(f, D, delimiter = ' ')
        
xgv = np.linspace(-20, 10, 900)
zgv = np.linspace(-20, 10, 901)

X, Z = np.meshgrid(xgv, zgv)
#y2 = twoscale(xgv, l1 = 3.1, l2 = 7.2)
y2 = threescale(xgv, l0 = 0.9, l1 = 3.1, l2 = 0.1)
D = np.zeros(X.shape)
for i in range(len(zgv)):
    D[i,:] = y2

fn = r"C:\Users\Scott\Documents\LSP\Water columns\Dual exponential tests\mycol2D.dat"
write2DLSP(fn, xgv, zgv, D)


plt.figure(1)
plt.clf()
plt.plot(xgv1, edens2[edens2.shape[0]/2]/2.65)
plt.plot(xgv1, edens1[edens1.shape[0]/2]/2.65)
plt.plot(xgv, y2)
plt.ylim(0.1e20, 1e24)
plt.hlines(1.71e21/2.65, -20, 10, linestyle='--')
#plt.ylim(0.01, 4)
#plt.hlines(nc, -20, 20)
plt.yscale('log')