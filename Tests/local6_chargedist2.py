# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:57:20 2016

@author: Scott
"""

import sys
try:
    import lspreader2 as rd
except:
    print "Modifying path to include LSPreader"
    readerpath = r'C:\Users\Scott\Documents\Programming\Python\lspreader'
    sys.path.append(readerpath) # Add the LSP reader as taking precedence after all but the current directory
    import lspreader2 as rd

import scottplots as sp
import lstools as ls
import AnalyzeAll as aa
import matplotlib.pyplot as plt
import numpy as np
import sftools as sf
from pext import pextanalysis as pa
import scipy.constants as sc
from pext import quantities
from scipy.signal import savgol_filter


def rebin(a, new_shape):
    """
    Copied directly from: http://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
    and https://gist.github.com/zonca/1348792
    ONLY WORKS WHEN INDICES ARE EVEN

    Resizes a 2d array by averaging or repeating elements, 
    new dimensions must be integral factors of original dimensions
    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged, 
        if the new shape is bigger array elements are repeated
    See Also
    --------
    resize : Return a new array with the specified shape.
    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])
    """
    M, N = a.shape
    print M, N
    m, n = new_shape
    print m, n
    if m<M:
        print m, M/m, n, N/n
        return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)

    
## CHARGE DENSITY CALCULATION
fns = [r"C:\Users\Scott\Documents\temp\oct2016\curtest_rho1\sclr1.p4", r"C:\Users\Scott\Documents\temp\oct2016\curtest_rho1\sclr820.p4"]
#fns = [r"C:\Users\Scott\Documents\temp\7-5 videos\1p5 test\sclr1007.p4"]
data = ls.scalars2D(fns, fld_ids = ['Rho', 'RhoN1', 'RhoN2', 'RhoN3',  'RhoN4', 'RhoN5', 'RhoN6', 'RhoN7', 'RhoN8', 'RhoN9', 'RhoN10',  'RhoN11'])

xgv = data['xgv']*1e4 # x values in microns
zgv = data['zgv']*1e4
dx = np.mean(np.diff(xgv))# dx in microns
dz = np.mean(np.diff(zgv))

# Net charge (in each cell)
data['RhoPos'] = 1 * data['RhoN11'] + 1 * data['RhoN2'] + 2 * data['RhoN3'] + 3 * data['RhoN4'] + 4 * data['RhoN5'] + 5 * data['RhoN6'] + 6 * data['RhoN7'] + 7 * data['RhoN8'] + 8 * data['RhoN9']
data['RhoNeg'] = 1 * data['RhoN10']
data['RhoNet'] = data['RhoPos'] - data['RhoNeg']

# Add up electrons, protons, and each species of ions (following list in the .lsp file for species numbers)
#pnet = np.mean(1 * data['RhoN11'] + 1 * data['RhoN2'] + 2 * data['RhoN3'] + 3 * data['RhoN4'] + 4 * data['RhoN5'] + 5 * data['RhoN6'] + 6 * data['RhoN7'] + 7 * data['RhoN8'] + 8 * data['RhoN9'], 0)
#nnet = np.mean(1 * data['RhoN10'], 0)

#qnet = pnet - nnet
#qnet = data['RhoNet'][1] - data['RhoNet'][0]
qnet = data['RhoNet'][1]
qnet2 = data['Rho'][1]*1e-6/1.6e-19

edens = data['RhoNeg'][1]
pdens = data['RhoPos'][1] # Positive charge density, in 

plt.figure(10)
plt.clf()
#plt.plot(xgv, qnet[qnet.shape[0]/2])
plt.plot(xgv, np.mean(qnet,0), alpha=0.2)
plt.plot(xgv, savgol_filter(np.mean(qnet,0), 81, 3))
plt.plot(xgv, np.mean(qnet2,0), alpha=0.2)

#plt.plot(xgv, nnet[qnet.shape[0]/2])
plt.ylim([-1e19, 1e19])
plt.xlim(-30, 10)
plt.axhline(y=0,linestyle="--", color='black')
#plt.yscale('log')

fig = sp.mypcolor(qnet, xgv, zgv, cmin=-0.2e21, cmax=0.2e21, cmap='RdBu', title='Qnet')
fig = sp.mypcolor(edens, xgv, zgv, cmax=1e21, title='Edens')
fig = sp.mypcolor(qnet2, xgv, zgv, cmin=-0.2e21, cmax=0.2e21, cmap='RdBu', title='Qnet2')

qnetsmall = rebin(qnet[0:400,0:400], (40,40))
sp.mypcolor(qnetsmall, xgv, zgv, cmin=-1e20, cmax=1e20, cmap='RdBu', title='Qnet Small')

qnet2small = rebin(qnet2[0:400,0:400], (40,40))
sp.mypcolor(qnet2small, xgv, zgv, cmin=-5e20, cmax=5e20, cmap='RdBu', title='Qnet2 Small')

#edenssmall = rebin(edens[0:960,0:960], (40,40))
#sp.mypcolor(edenssmall, xgv, zgv, cmax=1e21, title='Edens Small')