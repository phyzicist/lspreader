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

p4dir = r'C:\Users\Scott\Documents\temp\newcoltests\TNSA test2'
fns_fld, fns_scl = ls.getFldScl(p4dir)
data = ls.fields2D(fns_fld, fld_ids=['Ex', 'Ey', 'Ez', 'Vave10x', 'Vave10y', 'Vave10z','Vave11x', 'Vave11y', 'Vave11z'])
#data_scl = ls.scalars2D(fns_scl, fld_ids=['Vave11x', 'Vave11y', 'Vave11z'])


massE = 0.511 # electron rest mass, in MeV
massp = 938 # proton rest mass, in MeV

venorm = np.sqrt(data['Vave10x']**2 + data['Vave10y']**2 + data['Vave10z']**2)
#venorm = np.sqrt(data['Vave10x']**2 + data['Vave10y']**2 + data['Vave10z']**2)
elecKE =(np.sqrt(venorm**2 + 1) - 1)*massE # Electron kinetic energy, in MeV (kinetic energy of each constituent real electron, that is; not of the whole macroparticle)

vpnorm = np.sqrt(data['Vave11x']**2 + data['Vave11y']**2 + data['Vave11z']**2)
protKE =(np.sqrt(vpnorm**2 + 1) - 1)*massp # Electron kinetic energy, in MeV (kinetic energy of each constituent real electron, that is; not of the whole macroparticle)

KEeav = np.mean(elecKE,0)
KEpav = np.mean(protKE,0)
sp.mypcolor(KEeav, data['xgv'], data['zgv'])
sp.mypcolor(KEpav, data['xgv'], data['zgv'])

plt.figure()
plt.clf()
#plt.plot(data['xgv'], np.mean(KEeav[280:320,:],0))
plt.plot(data['xgv'], np.mean(KEpav[280:320,:],0))

#### DENSITY VISUAL TEST
#p4dir = r'C:\Users\Scott\Documents\temp\newcoltests\TNSA test2'
#data = ls.readFldScl(p4dir)
#sp.plotDens(data, outdir=p4dir, shortname = '', alltime=False)
#
#xgv = data['xgv']*1e4 # x values in microns
#zgv = data['zgv']*1e4
#dx = np.mean(np.diff(xgv))# dx in microns
#dz = np.mean(np.diff(zgv))
#
### CALCULATIONS
## Mean electron density
#edens = np.mean(data['RhoN10'],0)
#
#plt.figure(10)
##plt.clf()
#plt.plot(xgv, edens[edens.shape[0]/2])
#plt.yscale('log')
# PEXTTEXT
#p4dir = r'C:\Users\Scott\Documents\temp\mar1test\hres_osc'
#outdir = r'C:\Users\Scott\Documents\temp\mar1test\hres_osc'
#shortname = 'test'
#pextarr = pa.pextFull(p4dir, outdir = outdir, shortname = shortname, Utot_Jcm = 25.585283)


### PMOVIE TEST
#fn = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc\pmovie1.p4.gz"
#frames = rd.read_movie2(fn)
