# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 15:16:05 2016

@author: Scott
"""
import lspreader2 as rd

fn = r'C:\Users\Scott\Documents\temp\Jan 14 tests\5um_Foc-14 comparisons\Lres\sclr550.p4'
doms, header = rd.read_flds2(fn)

fld_ids = list(set(doms[0].keys()) - set(['xgv','ygv','zgv']))
flds = {}
for fld_id in fld_ids:
    flds[fld_id], xgv, zgv = rd.stitch2D(doms, fld_id)

import matplotlib as mpl

xgv = xgv * 1e4
zgv = zgv * 1e4

title = "Electron density, t = 110fs"
fignum = 2
C = flds['RhoN10']
fig = plt.figure(fignum)
plt.clf() # Clear the figure
cmin = np.min(C)
cmax = 1.74e21
ax = plt.subplot(111)
xr = [xgv[0],xgv[-1]] # min and max of xgv, needed for pcolorfast
zr = [zgv[0],zgv[-1]]
im = ax.pcolorfast(xr,zr,C, cmap='viridis')
ax.set_xlabel(r'X ($\mu m$)')
ax.set_ylabel(r'Z ($\mu m$)')
ax.set_title(title, fontsize=20)
cbar = fig.colorbar(im, label=r'Electron density (number/cm$^3$)')
im.set_clim(vmin=cmin, vmax=cmax)
plt.suptitle('Yellow is critical. 5 um scale length, $\lambda$/16 sim, focal point at X = +10um')
plt.axis('equal')
ax.set_xlim([-35,5])
ax.set_ylim([-20,20])

title = "Z=0 lineout"
fignum = 3
lineout = C[C.shape[0]/2,:]
fig = plt.figure(fignum)
ax = plt.subplot(111)
im = ax.plot(xgv,lineout)
ax.set_xlabel(r'X ($\mu m$)')
ax.set_ylabel(r'Density (number/cc)')
ax.set_title(title, fontsize=20)