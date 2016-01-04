# -*- coding: utf-8 -*-
"""
Open file, then do 2D of FFT on the file.

Created on Wed Dec 30 19:32:07 2015

@author: Scott
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np

fn = r'C:\Users\Scott\Documents\temp\lastest\run2_lastest longer\fields2D.hdf5'
f = h5py.File(fn, 'r')
times = f['times'][...]
fns = f['filenames'][...]
xgv = f['xgv'][...]*1e4
zgv = f['zgv'][...]*1e4
Ez = f['Ez'][...]

ix = np.argsort(times)

dx = xgv[1] - xgv[0]
dz = zgv[1] - zgv[0]
dt = times[ix[2]] - times[ix[1]]

# Some basic calculations
c = 3e8 # Speed of light in m/s
wl = 0.8e-6 # Wavelength of laser in m
fr_fund = c/wl # Frequency of the laser (the 'fundamental'), in Hz

freqHz = np.fft.rfftfreq(len(times), d = dt/1e9) # Frequencies in Hz, for upcoming real FFT
freq = freqHz/fr_fund # Frequencies in units of the fundamental

Imap = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map
FTmap1 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
FTmap2 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
pwr_sum = np.zeros(freq.shape)

map1_condit = np.logical_and(freq > 0.9, freq < 1.1)
map2_condit = np.logical_and(freq > 0.3, freq < 0.7)

for i in range(len(zgv)):
    for j in range(len(xgv)):
        tline = Ez[:,i,j]
        
        # Make a simple intensity map        
        esum = np.sum(tline**2)
        Imap[i,j] = esum
        
        # Make a more complex ftransform map
        time_ft = np.fft.rfft(tline[ix])
        pwr = np.abs(time_ft)**2
        FTmap1[i,j] = np.mean(pwr[map1_condit])
        FTmap2[i,j] = np.mean(pwr[map2_condit])
        pwr_sum = pwr_sum + pwr;
f.close()


# Make some plots
# Time-integrated intensity map (a.u.)
plt.figure(1)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(xgv,zgv,Imap, cmap='gray')
ax.set_xlabel('X (um)')
ax.set_ylabel('Z (um)')
ax.set_title('Integrated power from Ez (a.u.)')
X,Z = np.meshgrid(xgv,zgv)
plt.contour(X,Z,Imap)

# Fourier transform subset map 1
plt.figure(2)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(xgv,zgv,FTmap1, cmap='gray')
ax.set_xlabel('X (um)')
ax.set_ylabel('Z (um)')
ax.set_title('Fundamental map (0.9 to 1.1 x fundamental) (a.u.)')

# Fourier transform subset map 2
plt.figure(3)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(xgv,zgv,FTmap2, cmap='gray')
ax.set_xlabel('X (um)')
ax.set_ylabel('Z (um)')
ax.set_title('1/2-omega map (0.3 to 0.7 x fundamental) (a.u.)')


# Average power
plt.figure(4)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.plot(freq, np.log10(pwr_sum))
ax.set_xlabel('Frequency (normalized to laser fundamental)')
ax.set_ylabel('Log10(Power spectrum) (a.u.)')
ax.set_title('Power spectrum (a.u.)')
