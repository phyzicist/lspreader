# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 15:30:36 2015

@author: Scott
"""

# A test of the lspreader. Specifically, reading in the second

#import h5py
import lspreader2 as rd
import matplotlib.pyplot as plt
import numpy as np

fn = r"C:\Users\Scott\Documents\temp\lsp laser test fields\flds2090.p4"

doms2, header = rd.read_flds2(fn, flds=['E','B'])

fld, xgv, zgv = rd.stitch2D(doms2,'Ex')
fld_min, fld_max = -np.abs(fld).max(), np.abs(fld).max()

plt.figure(1)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(xgv*1e3,zgv*1e3,fld,cmap='RdBu', vmin=fld_min, vmax=fld_max)
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_title('Ex')
#plt.axis('equal')