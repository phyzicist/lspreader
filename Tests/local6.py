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

import scipy.constants as sc

def plotBackscat(data, trange = (80, 1.0e9)):
    # Assume equal spacing in time and space, and get the deltas
    dt = np.mean(np.diff(data['times']))*1e-9 # dt in seconds
    dx = np.mean(np.diff(data['xgv']))*1e-2 # dx in m
    dz = np.mean(np.diff(data['zgv']))*1e-2 # dx in m

    # Do the Fourier analysis
    maps, cuts, pwrsum, freq, df, nframes = sp.emFFT(data, kind='Backscat', trange=trange)
    
    # All these energies returned by emFFT are the mean energy value, across the time period; if we want the sum, we need to multiply by number of frames analyzed
    # But there is now another problem. We have overcounted the energy by a factor of (dx*c/dt), because while dx should be c/dt long, it is actually larger and a delta pulse of light gets counted multiple frames in a row.
    pwrline = pwrsum*(nframes * sc.c * dt/dx) # Units are now J per y meter per freq. (integrated rather than averaged over all time)
    lines = {}
    for k in maps:
        lines[k] = np.squeeze(maps[k]*(nframes * sc.c * dt/dx)) # Units are now Joules per y meter per z meter (integrated rather than averaged over all time). Not really, I can't figure out what's wrong though.
    
    # That's better!

    print lines['0'].shape
    
    
    print "Total energy, in mJ/micron:", np.sum(lines['all'])*dx*dz*1e-3
    print "Total energy, in mJ/micron:", np.sum(pwrline)*df*1e-3
    
    plt.figure(0)
    plt.plot(freq, pwrline*1e-3)
    plt.ylabel('mJ/ymicron/freq.')
    plt.xlabel('Freq.')

    plt.figure(1)
    plt.plot(data['zgv'], lines['all']*1e-9)
    plt.ylabel('mJ/yum/zum')
    plt.xlabel('Z ($\mu m$)')



    
#p4root = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump'
#outroot = r'C:\Users\Scott\Documents\temp\sclrtest\analysis'
#aa.analyzeAll(p4root, outroot)


#p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test0-2015-15-21_1923'
p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test2'
#data = ls.readFldScl(p4dir)
#
shortname = 'test0'
outdir = r'C:\Users\Scott\Documents\temp\sclrtest\analysis\test0'

plotBackscat(data, trange = (0, 60))
#chk_fs = 1.0#fs How many femtoseconds per (large-sized) data chunk for analysis?
#freqchunks = ls.chunkData(data, chk_fs) + ls.chunkData(data, chk_fs, offset_fs = chk_fs/2.0)
#denschunks = ls.chunkData(data, chk_fs/2.0) # Same number of chunks as bigchunks, but less long in time.
