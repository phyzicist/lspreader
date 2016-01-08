# -*- coding: utf-8 -*-
"""
Analyze a batch of p4 files for frequency components

Created on Wed Dec 30 15:09:37 2015

@author: Scott
"""
print "Importing functions"
import h5py
import os
import lspreader2 as rd
import numpy as np
from h5stitch2D import chunkIt, getfnsp4, fields2D, h5fields2D, h5fields2Dpar2, getTimes
import gzip
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def freqanalyze(fns):
	print "Will read ", len(fns), "files."
	# Load in the data
	divsp = 2 # Integer divisor by which to reduce the spatial resolution
	data = fields2D(fns, fld_ids = ['Ez'])
	print "Data read in."
	times = data['times']*1e6 # times, converted to fs
	fns = data['filenames']
	xgv = data['xgv'][::divsp]*1e4 # spatial X, converted to microns
	zgv = data['zgv'][::divsp]*1e4 # spatial Z, converted to microns
	Ez = data['Ez'][:,::divsp,::divsp]

	data2 = {}
	data2['times_fs'] = times
	data2['xgv_um'] = xgv
	data2['zgv_um'] = zgv
	data2['fns'] = fns

	print "times, fns, xgv, zgv, Ez"
	print times.shape
	print fns.shape
	print xgv.shape
	print zgv.shape
	print Ez.shape

	# Assume equal spacing in time and space, and get the deltas
	print "Calculating dt, dx, dz"
	dt = np.mean(np.diff(times))
	dx = np.mean(np.diff(xgv))
	dz = np.mean(np.diff(zgv))

	print dt, dx, dz

	# Calculate the frequency of the laser from its wavelength (800 nm)
	c = 3e8 # Speed of light in m/s
	wl = 0.8e-6 # Wavelength of laser in m
	fr_fund = c/wl # Frequency of the laser (the 'fundamental'), in Hz

	freqHz = np.fft.rfftfreq(len(times), d = dt/1e15) # Frequencies in Hz, for upcoming real FFT
	freq = freqHz/fr_fund # Frequencies in units of the fundamental

	data2['freq'] = freq

	# Make some interesting maps
	Imap = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map
	FTmap1 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
	FTmap2 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
	FTmap3 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
	FTmap4 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
	pwr_sum = np.zeros(freq.shape)

	map1_condit = np.logical_and(freq > 0.3, freq < 0.7)
	map2_condit = np.logical_and(freq > 0.9, freq < 1.1)
	map3_condit = np.logical_and(freq > 1.3, freq < 1.7)
	map4_condit = np.logical_and(freq > 1.8, freq < 2.3)

	print "Doing calculations"
	for i in range(len(zgv)):
		for j in range(len(xgv)):
			tline = Ez[:,i,j]
			
			# Make a simple intensity map        
			esum = np.sum(tline**2)
			Imap[i,j] = esum
			
			# Make a more complex ftransform map
			time_ft = np.fft.rfft(tline)
			pwr = np.abs(time_ft)**2
			FTmap1[i,j] = np.sum(pwr[map1_condit])
			FTmap2[i,j] = np.sum(pwr[map2_condit])
			FTmap3[i,j] = np.sum(pwr[map3_condit])
			FTmap4[i,j] = np.sum(pwr[map4_condit])
			pwr_sum = pwr_sum + pwr;



	data2['pwr_sum'] = pwr_sum
	data2['Imap'] = Imap
	data2['FTmap_0_5'] = FTmap1
	data2['FTmap_1'] = FTmap2
	data2['FTmap_1_5'] = FTmap3
	data2['FTmap_2'] = FTmap4

	return data2

def plotme(data2):
	xgv = data2['xgv_um']
	zgv = data2['zgv_um']
	freq = data2['freq']
	times = data2['times_fs']
	pwr_sum = data2['pwr_sum']
	Imap = data2['Imap']
	FTmap1 = data2['FTmap_0_5']
	FTmap2 = data2['FTmap_1']
	FTmap3 = data2['FTmap_1_5']
	FTmap4 = data2['FTmap_2']

	meantime = np.mean(times)

	# Make some plots
	print "Making plots."
	figs = [None]*6
	# Time-integrated intensity map (a.u.)
	figs[0] = plt.figure(1)
	plt.clf() # Clear the figure
	ax = plt.subplot(111)
	ax.pcolorfast(xgv,zgv,Imap, cmap='gray')
	ax.set_xlabel('X (um)')
	ax.set_ylabel('Z (um)')
	ax.set_title('Integrated power from Ez (a.u.)')
	X,Z = np.meshgrid(xgv,zgv)
	plt.contour(X,Z,Imap)
	ax.text(-19, 17, r't=' + "{:.1f}".format(meantime) + " fs", fontsize=24, color='white')

	# Average power
	figs[1] = plt.figure(2)
	plt.clf() # Clear the figure
	ax = plt.subplot(111)
	ax.plot(freq, np.log10(pwr_sum))
	ax.set_xlabel('Frequency (normalized to laser fundamental)')
	ax.set_ylabel('Log10(Power spectrum) of Ez (a.u.)')
	ax.set_title('Log of Power spectrum) at t=' + "{:.1f}".format(meantime) + " fs")

	# Fourier transform subset map 1
	figs[2] = plt.figure(3)
	plt.clf() # Clear the figure
	ax = plt.subplot(111)
	ax.pcolorfast(xgv,zgv,FTmap1, cmap='gray')
	ax.set_xlabel('X (um)')
	ax.set_ylabel('Z (um)')
	ax.set_title('Ez, Half omega map (0.3 to 0.7 x fundamental) (a.u.)')
	ax.text(-19, 17, r't=' + "{:.1f}".format(meantime) + " fs", fontsize=24, color='white')
	
	# Fourier transform subset map 2
	figs[3] = plt.figure(4)
	plt.clf() # Clear the figure
	ax = plt.subplot(111)
	ax.pcolorfast(xgv,zgv,FTmap2, cmap='gray')
	ax.set_xlabel('X (um)')
	ax.set_ylabel('Z (um)')
	ax.set_title('Ez, Fundamental map (0.9 to 1.1 x fundamental) (a.u.)')
	ax.text(-19, 17, r't=' + "{:.1f}".format(meantime) + " fs", fontsize=24, color='white')

	# Fourier transform subset map 3
	figs[4] = plt.figure(5)
	plt.clf() # Clear the figure
	ax = plt.subplot(111)
	ax.pcolorfast(xgv,zgv,FTmap3, cmap='gray')
	ax.set_xlabel('X (um)')
	ax.set_ylabel('Z (um)')
	ax.set_title('Ez, Three-halves omega map (1.3 to 1.7 x fundamental) (a.u.)')
	ax.text(-19, 17, r't=' + "{:.1f}".format(meantime) + " fs", fontsize=24, color='white')

	# Fourier transform subset map 4
	figs[5] = plt.figure(6)
	plt.clf() # Clear the figure
	ax = plt.subplot(111)
	ax.pcolorfast(xgv,zgv,FTmap4, cmap='gray')
	ax.set_xlabel('X (um)')
	ax.set_ylabel('Z (um)')
	ax.set_title('Ez, Two omega map (1.8 to 2.3 x fundamental) (a.u.)')
	ax.text(-19, 17, r't=' + "{:.1f}".format(meantime) + " fs", fontsize=24, color='white')

	print "Saving figures"
	for i in range(len(figs)):
		figs[i].savefig('fig' + str(i) + '.png')

	return figs

def h5data2(data2, fn = 'data2.hdf5'):
	"""Save the data2 dictionary to an HDF5 file, assuming it contains only NumPy arrays"""
	with h5py.File(fn, 'w') as f:
		for k in data2:
			f.create_dataset(k, data = data2[k], compression='gzip', compression_opts=4)

if __name__ == "__main__":
	## MAIN PROGRAM
	#folder = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir'
	folder = r'/tmp/ngirmang.1-2d-nosolid-151113/'
	print "Getting files in folder"
	fns = getfnsp4(folder)[400:440]
	#h5fields2D(folder, fld_ids = 'Ez'])
	#h5fields2Dpar2(folder, fld_ids = ['Ez'])
	data2 = freqanalyze(fns)
	print "Saving HDF5 file."
	figs = plotme(data2)
