# -*- coding: utf-8 -*-
"""
Open file, then do 2D of FFT on the file.

Created on Wed Dec 30 19:32:07 2015

@author: Scott
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack

fn = r'C:\Users\Scott\Documents\temp\lastest\fields2D.hdf5'
f = h5py.File(fn, 'r')
times = f['times'][...]
fns = f['filenames'][...]
xgv = f['xgv'][...]*1e4
zgv = f['zgv'][...]*1e4
Ex = f['Ez'][...]

ix = np.argsort(times)

dx = xgv[1] - xgv[0]
dz = zgv[1] - zgv[0]
dt = times[ix[10]] - times[ix[9]]

Imap = np.zeros((Ex.shape[1],Ex.shape[2]))
for i in range(len(zgv)):
    for j in range(len(xgv)):
        tline = Ex[:,i,j]
        esum = np.sum(tline**2)
        Imap[i,j] = esum

fld = f['Ez'][ix[91]]

f.close()


# Make some plots
plt.figure(1)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(xgv,zgv,Imap, cmap='gray')
ax.set_xlabel('X (um)')
ax.set_ylabel('Z (um)')
ax.set_title('Integrated power from Ez (a.u.)')
