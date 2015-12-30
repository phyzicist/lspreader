# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 15:30:36 2015

@author: Scott
"""

# A test of the lspreader

#import h5py
import lspreader2 as rd
import matplotlib.pyplot as plt
import numpy as np

fn = r"C:\Users\Scott\Documents\temp\lsp laser test fields\flds2090.p4"

data, header, doms2 = rd.read_flds2(fn, flds=['E','B'])

z = np.squeeze(doms2[0]['Ex'])
zX = np.squeeze(doms2[0]['x'])
zY = np.squeeze(doms2[0]['y'])
zZ = np.squeeze(doms2[0]['z'])

for i in range(1,48):
    z_tmp = np.squeeze(doms2[i]['Ex'])[1:,:]
    zX_tmp = np.squeeze(doms2[i]['x'])[1:,:]
    zY_tmp = np.squeeze(doms2[i]['y'])[1:,:]
    zZ_tmp = np.squeeze(doms2[i]['z'])[1:,:]
    z = np.concatenate((z,z_tmp),0)
    zX = np.concatenate((zX,zX_tmp),0)
    zY = np.concatenate((zY,zY_tmp),0)
    zZ = np.concatenate((zZ,zZ_tmp),0)

z_min, z_max = -np.abs(z).max(), np.abs(z).max()

zgv = zZ[:,0]
xgv = zX[0,:]
x,y = np.meshgrid(xgv,zgv)

plt.figure(1)
plt.clf() # Clear the figure
ax = plt.subplot(111)
ax.pcolorfast(xgv*1e3,zgv*1e3,z,cmap='RdBu', vmin=z_min, vmax=z_max)
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_title('Ex')
#plt.axis('equal')

plt.figure(2)
plt.clf() # Clear the figure
ax = plt.subplot(4, 1, 1)
#ax.pcolorfast(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
ax.pcolorfast(z, cmap='RdBu', vmin=z_min, vmax=z_max)

ax2 = plt.subplot(4, 1, 2)
#ax.pcolorfast(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
ax2.pcolorfast(zX, cmap='RdBu', vmin=np.min(zX), vmax=np.max(zX))

ax2 = plt.subplot(4, 1, 3)
#ax.pcolorfast(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
ax2.pcolorfast(zZ, cmap='RdBu', vmin=np.min(zZ), vmax=np.max(zZ))