# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 15:51:24 2016

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


def getstuff(fn):
    """ Get ionization data """
    data = ls.scalars2D([fn])
    
    xgv = data['xgv']*1e4 # x values in microns
    zgv = data['zgv']*1e4
    dx = np.mean(np.diff(xgv))# dx in microns
    dz = np.mean(np.diff(zgv))
    
    ## CALCULATIONS
    # Mean electron density
    edens = np.mean(data['RhoN10'],0)
    
    # Mean ion density
    pdens = np.mean(data['RhoN11'],0)
    
    # Mean oxygen ionization state
    old_settings = np.seterr(divide='ignore', invalid='ignore') # Set to ignore divide by zero error. (We will divide by zero where no ions exist)
    ionstate = np.mean((0*data['RhoN1'] + 1*data['RhoN2'] + 2*data['RhoN3'] + 3*data['RhoN4'] + 4*data['RhoN5'] + 5*data['RhoN6'] + 6*data['RhoN7'] + 7*data['RhoN8'])/(data['RhoN1'] + data['RhoN2'] + data['RhoN3'] + data['RhoN4'] + data['RhoN5'] + data['RhoN6'] + data['RhoN7'] + data['RhoN8']), 0)
    ionstate = np.nan_to_num(ionstate)
    np.seterr(**old_settings)  # reset divide by zero error to prior settings
    
    return xgv, zgv, ionstate, edens, pdens, data

fn1 = r'C:\Users\Scott\Documents\temp\mar1test\ioninvest\sclr1.p4.gz'
fn2 = r'C:\Users\Scott\Documents\temp\mar1test\ioninvest\sclr1734.p4.gz'

xgv, zgv, ionstate1, edens1, pdens1, data1 = getstuff(fn1)
_, _, ionstate2, edens2, pdens2, data2 = getstuff(fn2)

X, Z = np.meshgrid(xgv, zgv)

#fig = plt.figure(1)
#plt.clf()
#plt.title("Ionization state")
#plt.pcolor(X, Z, ionstate2, vmin=5, vmax=7)
#plt.colorbar()
#
#fig = plt.figure(2)
#plt.clf()
#plt.title("Electron density 2")
#plt.pcolor(X, Z, edens1, vmin = 0, vmax = 2.0e21)
#plt.colorbar()
#
#fig = plt.figure(3)
#plt.clf()
#plt.title("Electron density 2")
#plt.pcolor(X, Z, edens2, vmin = 0, vmax = 2.0e21)
#plt.colorbar()

eratio = edens2/edens1
fig = plt.figure(4)
plt.clf()
plt.title("Electron density final per initial")
plt.pcolor(X, Z, eratio, vmin = 2.65, vmax = 2.65+0.33, cmap='magma')
plt.colorbar()