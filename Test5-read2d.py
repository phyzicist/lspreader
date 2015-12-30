# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 16:51:40 2015

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

ix = np.argsort(times)

fld = f['Ez'][ix[91]]
fld_min, fld_max = -np.abs(fld).max(), np.abs(fld).max()

plt.figure(1)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(xgv*1e4,zgv*1e4,fld,cmap='RdBu', vmin=fld_min, vmax=fld_max)
ax.set_xlabel('X (um)')
ax.set_ylabel('Z (um)')
ax.set_title('Ex' + str(times[ix[91]]*1e6))
#plt.axis('equal')

plt.figure(2)
plt.clf()
plt.plot(times[ix])

#plt.figure(3)
#plt.clf()

#ax = plt.subplot(111)
#ax.pcolorfast(np.fft.fft2(fld))
dx = xgv[1] - xgv[0]
dz = zgv[1] - zgv[0]
dt = times[ix[10]] - times[ix[9]]

# Take the fourier transform of the image.
F1 = np.fft.rfft2(fld)
F2 = np.fft.fftshift(F1)
# Calculate a 2D power spectrum
psd2D = np.log10(np.abs(F2)**2)
freqz = np.fft.rfftfreq(fld.shape[0], d = dz)
freqx = np.fft.rfftfreq(fld.shape[1], d = dx)
plt.figure(3)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(freqx[0:150], freqz[0:150], psd2D[0:150,0:190])
ax.set_title('Ex 2D FFT')
#plt.axis('equal')


lnt = fld[200,:]
tr1D = np.fft.rfft(lnt)

pwr = np.log10(np.abs(tr1D)**2)

plt.figure(4)
plt.clf()
ax = plt.subplot(111)
#xvals = np.array(range(len(tr1D)))*dx
freq = np.fft.rfftfreq(len(lnt), d = dx)
idx = np.argsort(freq)
ax.plot(freq[idx],pwr[idx])
ax.set_title('Power spectrum for lineout')
ax.set_xlabel('Spatial frequency (cycles/micron)')

plt.figure(5)
plt.clf()
ax = plt.subplot(111)
ax.plot(xgv,lnt)
ax.set_title('Ex, central lineout')

f.close()