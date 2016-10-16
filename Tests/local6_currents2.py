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

import os
import scottplots as sp
import lstools as ls
import AnalyzeAll as aa
import matplotlib.pyplot as plt
import numpy as np
import sftools as sf
from pext import pextanalysis as pa
import scipy.constants as sc
from pext import quantities
from mpl_toolkits.mplot3d import Axes3D

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

def readFldScl2(p4dir, fld_ids=['Ex','Ey','Ez','Bx','By','Bz','Jx','Jy','Jz'], scl_ids=['Rho', 'RhoN1', 'RhoN2', 'RhoN3',  'RhoN4', 'RhoN5', 'RhoN6', 'RhoN7', 'RhoN8', 'RhoN9', 'RhoN10',  'RhoN11'], divsp=1, divt=1, pool=None):
    """ Read matching field and scalar files (default fld_ids) in a directory into a single data array """
    ## READ MATCHING FIELDS AND SCALARS INTO DATA ARRAY
    fns_fld, fns_scl = ls.getFldScl(p4dir)
    
    data_fld = ls.fields2D(fns_fld[::divt], fld_ids=fld_ids, divsp=divsp, pool=pool)
    data_scl = ls.scalars2D(fns_scl[::divt], fld_ids=scl_ids, divsp=divsp, pool=pool)
    
    # Sanity check that times, xgv, and zgv are essentially identical between the two
    if np.max(np.abs(data_scl['xgv'] - data_fld['xgv'])) + np.max(np.abs(data_scl['zgv'] - data_fld['zgv'])) + np.max(np.abs(data_scl['times'] - data_fld['times'])) > 0.000001:
        raise Exception("Scalar and field files grid vectors (or times) do not match!! This shouldn't be the case. Major problem with file choice or with lspreader.")
    
    data = {} # The merged data array. Note that no copies of arrays will be made, just pointers shuffled around.
    for k in data_fld.keys():
        data[k] = data_fld[k]
    for k in set(data_scl.keys()).difference(data_fld.keys()): # Add only those keys not already in the key list
        data[k] = data_scl[k]
    return data
#p4root = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump'
#outroot = r'C:\Users\Scott\Documents\temp\sclrtest\analysis'
#aa.analyzeSome(p4root, outroot, ['a50f-14_mres_so'])


#p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test0-2015-15-21_1923'
#p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test2'
#data = ls.readFldScl(p4dir)
#
#shortname = 'test0'
#outdir = r'C:\Users\Scott\Documents\temp\sclrtest\analysis\test0'
#pa.pextFull(p4dir, outdir=outdir, shortname=shortname)
#sp.plotBScat(data, tcut=65, outdir=outdir, shortname = shortname)


#chk_fs = 1.0#fs How many femtoseconds per (large-sized) data chunk for analysis?
#freqchunks = ls.chunkData(data, chk_fs) + ls.chunkData(data, chk_fs, offset_fs = chk_fs/2.0)
#denschunks = ls.chunkData(data, chk_fs/2.0) # Same number of chunks as bigchunks, but less long in time.

full=False

shortname = r'bext_slabcomb'
#shortname = r'curtest_240fs_slabcomb'
p4dir = os.path.join(r'C:\Users\Scott\Documents\temp\oct2016', shortname)
if full:
    data = readFldScl2(p4dir)

    xgv = data['xgv']*1e4 # x values in microns
    zgv = data['zgv']*1e4
    dx = np.mean(np.diff(xgv))# dx in microns
    dz = np.mean(np.diff(zgv))


i=27

tstring = 't=' + "{:.1f}".format(data['times'][i]*1e6)

edens = data['RhoN10'][i]
pdens = data['RhoN11'][i]
Jzdens = data['Jz'][i]
Jxdens = data['Jx'][i]

fig = plt.figure(7)
C = edens
(cmin, cmax) = (0, 3e21)
sticker = '$e^-$'
title = 'Electron density'
clabel = 'Density (number/cc)'
fld_id = r'$\rho$'
sp.mypcolor(C, xgv, zgv, cmin=cmin,  cmax=cmax, fig=fig, title=title, tstring=tstring, clabel=clabel, fld_id=fld_id, sticker=sticker, rfooter=shortname, edens=edens)

fig = plt.figure(20)
C = pdens
(cmin, cmax) = (0, 3e20)
sticker = '$p^+$'
title = 'Proton density'
clabel = 'Density (number/cc)'
fld_id = r'$\rho$'
sp.mypcolor(C, xgv, zgv, cmin=cmin,  cmax=cmax, fig=fig, title=title, tstring=tstring, clabel=clabel, fld_id=fld_id, sticker=sticker, rfooter=shortname)


fig = plt.figure(1)
C = Jxdens
(cmin, cmax) = (-8e11, 8e11)
cmap = 'RdBu'
sticker = 'all'
title = 'Current density, Jx'
clabel = 'Current density (amps/cm^2)'
fld_id = r'$Jx$'
fig = sp.mypcolor(C, xgv, zgv, cmin=cmin,  cmax=cmax, cmap=cmap, fig=fig, title=title, tstring=tstring, clabel=clabel, fld_id=fld_id, sticker=sticker, rfooter=shortname)

fig = plt.figure(2)
fig = sp.mypcolor(Jzdens, xgv, zgv, fig=fig, cmin=-8e11, cmax=8e11, cmap='RdBu', title='Jzdens')

fig = plt.figure(3)
Jzdenssmall = rebin(Jzdens[0:140,0:400], (20,40))
sp.mypcolor(Jzdenssmall, xgv, zgv, fig=fig, cmin=-8e11, cmax=8e11, cmap='RdBu', title='Jxdens Small')


Ez = data['Ez'][i]
fig = plt.figure(4)
fig = sp.mypcolor(Ez, xgv, zgv, fig=fig, cmin=-3e7, cmax=3e7, cmap='RdBu', title="Ez")

Ex = data['Ex'][i]
fig = plt.figure(21)
fig = sp.mypcolor(Ex, xgv, zgv, fig=fig, cmin=-1e7, cmax=1e7, cmap='RdBu', title="Ex")

By = data['By'][i]
fig = plt.figure(5)
fig = sp.mypcolor(By, xgv, zgv, cmin=-1e8, cmax=1e8, fig=fig, cmap='RdBu', title="By")

rho = data['Rho'][i] # Charge density, from microcoulombs/cm^3
rhoSI = rho * 1e-6 * (1e2)**3 # Charge density, in Coulombs/m^3
rhoSIsmall = rebin(rhoSI[0:140,0:400], (70,100))
xgvsmall = np.linspace(xgv[0], xgv[399], 100)
zgvsmall = np.linspace(zgv[0], zgv[139], 70)
Xs, Zs = np.meshgrid(xgvsmall, zgvsmall)

Vs = getPot(Xs*1e-6, Zs*1e-6, Xs*1e-6, Zs*1e-6, rhoSIsmall, xref=+5.0e-6)
fig = plt.figure(6)
sp.mypcolor(Vs/1e6, xgvsmall, zgvsmall, fig=fig, cmin=-0.5, cmax=5.5, title="Voltage (MV)", cmap='viridis')
plt.axis("equal")

fig = plt.figure(8)
plt.clf()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(Xs, Zs, -Vs/1e6, cmap='viridis')
ax.set_zlim(-5.0, 0.2)
ax.set_zlabel("Well depth (MeV)")
ax.set_xlabel("X (um)")
ax.set_ylabel("Z (um)")
ax.set_title("Electron potential well")

fig = plt.figure(9)
fig = sp.mypcolor(rhoSI, xgv, zgv, cmin=-2e7, cmax=2e7, fig=fig, cmap='RdBu', title='Rho')

# PEXTTEXT
#p4dir = r'C:\Users\Scott\Documents\temp\mar1test\hres_osc'
#outdir = r'C:\Users\Scott\Documents\temp\mar1test\hres_osc'
#shortname = 'test'
#pextarr = pa.pextFull(p4dir, outdir = p4dir, shortname = '', Utot_Jcm = 25.585283)


### PMOVIE TEST
#fn = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc\pmovie1.p4.gz"
#frames = rd.read_movie2(fn)
