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
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import cumtrapz


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

def myfunc(xval, yval, X, Y, R0, L):
    Ri = np.sqrt((X - xval)**2 + (Y - yval)**2) # Radius with respect to this point
    ct = (Ri > 0)
    return np.sum((L[ct]/(2*np.pi*sc.epsilon_0)) * np.log(R0[ct]/Ri[ct])) # Potential at this particular point

def getPot(Xv, Yv, X, Y, Rho, xref=-30.0e-6, yref=0.0):
    """ Get the potential at Xv, Yv points (all inputs and outputs in SI units) """
    # Convert to linear charge density pixelmap
    dA = np.mean(np.diff(X[0,:]))*np.mean(np.diff(Y[:,0])) # Assuming uniform grid, the area of one X/Y/Rho cell (m^2)
    print dA
    L = Rho * dA # Linear charge density (Coulombs/m) for each X/Y/Rho cell

    R0 = np.sqrt((X - xref)**2 + (Y - yref)**2) # Radius of each pixel with respect to the designated reference point

    # Numerically solve
    vfunc = np.vectorize(myfunc, excluded=[2,3,4,5])
    Vv = vfunc(Xv, Yv, X, Y, R0, L)

    return Vv

## CHARGE DENSITY CALCULATION
#fns = [r"C:\Users\Scott\Documents\temp\oct2016\curtest_7p2\sclr1.p4.gz", r"C:\Users\Scott\Documents\temp\oct2016\curtest_7p2\sclr989.p4.gz"]
fns = [r"C:\Users\Scott\Documents\temp\oct2016\curtest_7p2\sclr989.p4.gz"]
data = ls.scalars2D(fns, fld_ids = ['Rho', 'RhoN1', 'RhoN2', 'RhoN3',  'RhoN4', 'RhoN5', 'RhoN6', 'RhoN7', 'RhoN8', 'RhoN9', 'RhoN10',  'RhoN11'])
fns2 = [r"C:\Users\Scott\Documents\temp\oct2016\curtest_7p2\flds989.p4.gz"]
data2 = ls.fields2D(fns2)

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
qnet = data['RhoNet'][0]
qnet2 = data['Rho'][0]*1e-6/1.6e-19
rho = data['Rho'][0] # Charge density, from microcoulombs/cm^3
rhoSI = rho * 1e-6 * (1e2)**3 # Charge density, in Coulombs/m^3

edens = data['RhoNeg'][0]
pdens = data['RhoPos'][0] # Positive charge density, in 

Ex = data2['Ex'][0]
Ez = data2['Ez'][0]
Xt = -cumtrapz(Ex, xgv)
Zt = -cumtrapz(Ez[:,0], zgv)
Vsynth = Xt
Vsynth = (Xt[1:,:].T + Zt).T
Vsynth = Vsynth - np.min(Vsynth)
sp.mypcolor(Vsynth, xgv[1:], zgv[1:], cmin=0, cmax=np.max(np.abs(Vsynth)), cmap='viridis', title='VSynth')

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

rhoSIsmall = rebin(rhoSI[0:400,0:400], (100,100))
xgvsmall = np.linspace(xgv[0], xgv[399], 100)
zgvsmall = np.linspace(zgv[0], zgv[399], 100)

sp.mypcolor(qnetsmall, xgvsmall, zgvsmall, cmin=-1e20, cmax=1e20, cmap='RdBu', title='Qnet Small')

sp.mypcolor(rhoSIsmall, xgvsmall, zgvsmall, cmin=-5e20, cmax=5e20, cmap='RdBu', title='Rho (SI) Small')

Xs, Zs = np.meshgrid(xgvsmall, zgvsmall)
X, Z = np.meshgrid(xgv, zgv)
Vline = getPot(xgv*1e-6, np.min(np.abs(zgv))*1e-6, X*1e-6, Z*1e-6, rhoSI, xref=-40.0e-6)

Vs = getPot(Xs*1e-6, Zs*1e-6, Xs*1e-6, Zs*1e-6, rhoSIsmall, xref=-40.0e-6)
sp.mypcolor(Vs/1e6, xgvsmall, zgvsmall, title="Voltage (MV)")
#edenssmall = rebin(edens[0:960,0:960], (40,40))
#sp.mypcolor(edenssmall, xgv, zgv, cmax=1e21, title='Edens Small')


fig = plt.figure(2)
plt.clf()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(Xs, Zs, -Vs/1e6, cmap='viridis')
ax.set_zlim(-2.0, 0.2)
ax.set_zlabel("Well depth (MeV)")
ax.set_xlabel("X (um)")
ax.set_ylabel("Z (um)")
ax.set_title("Electron potential well")

fig = plt.figure(3)
plt.clf()
plt.plot(xgv, -Vline*1e-6)