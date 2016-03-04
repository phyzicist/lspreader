# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 17:21:46 2016

Version 4 attempts a read of a folder in serial

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
import matplotlib as mpl
from pmov import traj2

## PMOVIE TEST
fn1 = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc\pmovie1.p4.gz"
frames1 = rd.read_movie2(fn1)
dat1 = frames1[0]['data']

dat1s, stats1 = traj2.sortOne(fn1)

fn2 = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc\pmovie1250.p4.gz"
dat2s, stats2 = traj2.sortOne(fn2)

data_new = traj2.fillGaps(dat2s, dat1s)

# TODO: Save the dat3s slice into the HDF5 file
mydir = "C:\Users\Scott\Documents\temp\mar1test\hres_osc"
traj2.serTraj(mydir)

dath = dat2s

d = {}

massE = 0.511 # electron rest mass, in MeV
u_norm = np.sqrt(dath['ux']**2+dath['uy']**2+dath['uz']**2) # 'uz', etc. is actually momentum / (speed of light * rest mass); it's what they call gamma-beta units; that is, ux = px/c/mass_0 = gamma * beta_x
#p_norm = massE * u_norm # The momenta of the particles, 
d['KE'] =(np.sqrt(u_norm**2+1)-1)*massE # Electron kinetic energy, in MeV (kinetic energy of each constituent real electron, that is; not of the whole macroparticle)
d['q'] = -dath['q']*1e3 # Amount of negative charge, nanoColoumbs/cm
d['phi'] = np.arctan2(dath['uz'], dath['ux'])

pc_xdims = [-30, 10]
pc_zdims = [-15, 15]

ct =np.logical_and(d['KE'] > 0.200, np.abs(d['phi']) > np.deg2rad(180 - 40))
dath = dath[ct]# Part of data that will participate in histogramming
H, xedges, yedges = np.histogram2d(dath['x']*1e4, dath['z']*1e4, weights=-dath['q'], range = [pc_xdims, pc_zdims], bins=(40,50))

#Hx, xedges2 = np.histogram(dath['xi']*1e4, weights=-dath['q'], range=pc_xdims, bins=100)

fig = plt.figure(1)
plt.clf()
ax = fig.add_subplot(111)
ax.set_title('pcolormesh: exact bin edges')
X, Y = np.meshgrid(xedges, yedges)
im = ax.pcolormesh(X, Y, H.swapaxes(0,1), cmap='viridis')
plt.colorbar(im)
plt.axis('equal')

fig = plt.figure(2)
plt.clf()
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(dath['x']*1e4, 50, normed=1, range=pc_xdims, facecolor='green', alpha=0.75, weights=-dath['q'])