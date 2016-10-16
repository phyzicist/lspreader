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

if full:
    p4dir = r'C:\Users\Scott\Documents\temp\oct2016\curtest_grid2wire'
    fns_fld, fns_scl = ls.getFldScl(p4dir)
    data = ls.fields2D(fns_fld, fld_ids = ['Jx', 'Jy', 'Jz', 'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'])
    
    xgv = data['xgv']*1e4 # x values in microns
    zgv = data['zgv']*1e4
    dx = np.mean(np.diff(xgv))# dx in microns
    dz = np.mean(np.diff(zgv))


i=6
Jxdens = data['Jx'][i]
fig = plt.figure(1)
fig = sp.mypcolor(Jxdens, xgv, zgv, fig=fig, cmin=-8e10, cmax=8e10, cmap='RdBu', title='Jxdens')

Jzdens = data['Jz'][i]
fig = plt.figure(2)
fig = sp.mypcolor(Jzdens, xgv, zgv, fig=fig, cmin=-8e10, cmax=8e10, cmap='RdBu', title='Jzdens')

fig = plt.figure(3)
Jzdenssmall = rebin(Jzdens[0:400,0:400], (40,40))
sp.mypcolor(Jzdenssmall, xgv, zgv, fig=fig, cmin=-8e11, cmax=8e11, cmap='RdBu', title='Jxdens Small')


Ez = data['Ez'][i]
fig = plt.figure(4)
fig = sp.mypcolor(Ez, xgv, zgv, fig=fig, cmin=-3e7, cmax=3e7, cmap='RdBu', title="Ez")

By = data['By'][i]
fig = plt.figure(5)
fig = sp.mypcolor(By, xgv, zgv, cmin=-3e7, cmax=3e7, fig=fig, cmap='RdBu', title="By")

# PEXTTEXT
#p4dir = r'C:\Users\Scott\Documents\temp\mar1test\hres_osc'
#outdir = r'C:\Users\Scott\Documents\temp\mar1test\hres_osc'
#shortname = 'test'
#pextarr = pa.pextFull(p4dir, outdir = p4dir, shortname = '', Utot_Jcm = 25.585283)


### PMOVIE TEST
#fn = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc\pmovie1.p4.gz"
#frames = rd.read_movie2(fn)
