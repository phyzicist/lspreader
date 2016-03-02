# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 11:16:47 2016

Secondary analysis of realistic pre-pulse FLASH handoff sims (z batch).

@author: Scott
"""

import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

analydir = r'C:\Users\Scott\Box Sync\Analysis3'
namelist = next(os.walk(analydir))[1] # Looks for all folders in this list

goodnames = []
focs = []
for name in namelist:
    m1 = re.search('^zmar31f(-*[0-9]*?)$', name)
    if m1 and name != 'zmar31f5':
        focs.append(float(m1.groups()[0]))
        goodnames.append(name)

ix_sort = np.argsort(focs) # Ordering of indices that will go in order of focal depth
goodnames = np.array(goodnames)[ix_sort]
focs = np.array(focs)[ix_sort]

print goodnames
print focs

plas = np.ones(len(goodnames))
effics = []
Rlight = []
Rgreen = []
Rhalf = []
Romega = []
Rblue = []
effics2 = []
chargeall = []
chargehigh = []

for name in goodnames:
    fn_effic = os.path.join(analydir, name, name + ' - Electron Efficiency.csv')
    arr_effic = np.genfromtxt(fn_effic, skip_header=2, delimiter=',')
    effics.append(arr_effic[2,1])
    effics2.append(arr_effic[0,1] - arr_effic[2,1])
    chargehigh.append(arr_effic[3,3])
    chargeall.append(arr_effic[0,3])
    fn_refl = os.path.join(analydir, name, name + ' - Reflectivity.csv')
    try:
        arr_refl = np.genfromtxt(fn_refl, skip_header=1, delimiter=',')
        Rlight.append(arr_refl[0,1] - arr_refl[2,1]) # All reflectivity energy, less the static
        Rgreen.append(arr_refl[3,1])
        Rblue.append(arr_refl[1,1])
        Romega.append(arr_refl[5,1])
        Rhalf.append(arr_refl[6,1])
    except:
        Rlight.append(0)
        Rgreen.append(0)
        Rblue.append(0)
        Romega.append(0)
        Rhalf.append(0)



focs = np.array(focs)
effics = np.array(effics)
effics2 = np.array(effics2)
Rlight = np.array(Rlight)
Rgreen = np.array(Rgreen)
Rblue = np.array(Rblue)
Romega = np.array(Romega)
Rhalf = np.array(Rhalf)
chargeall = np.array(chargeall)
chargehigh = np.array(chargehigh)

print effics
print Rlight


# EFFICIENCY FIGURE
nrows = 2
ncols = 1
fig = plt.figure(1)
plt.clf()
ax = plt.subplot(nrows, ncols, 1)
plt.plot(focs, effics, label='High-energy efficiency')
plt.ylim(0,np.max(effics))
plt.title('High-energy conversion efficiency, $KE>500 keV, \pm 40^{\circ}$', fontsize=14)
plt.xlabel('Laser focus, distance past column surface ($\mu m$)', fontsize=14)

ax = plt.subplot(nrows, ncols, 2)
plt.plot(focs, effics2, label='Low-energy efficiency')
plt.ylim(0,np.max(effics2))
plt.title('Low-energy conversion efficiency, $120<KE<500 keV, \pm 40^{\circ}$', fontsize=14)
plt.xlabel('Laser focus, distance past column surface ($\mu m$)', fontsize=14)

plt.tight_layout()


# CHARGE FIGURE
nrows = 1
ncols = 1
fig = plt.figure(2)
plt.clf()

ax = plt.subplot(nrows, ncols, 1)
plt.plot(focs, chargeall, label='$KE>120 keV, \pm 40^{\circ}$')
#plt.plot(focs, chargemid, label='$KE>500 keV, \pm 40^{\circ}$')
plt.plot(focs, chargehigh, label='$KE>1 MeV, \pm 40^{\circ}$')
plt.ylim(0,np.max(chargeall))
plt.title('Backward charge', fontsize=14)
plt.ylabel('3d-equivalent charge, pC', fontsize=14)
plt.xlabel('Laser focus, distance past column surface ($\mu m$)', fontsize=14)
plt.legend()

plt.tight_layout()


# BACKSCATTER FIGURE
nrows = 3
ncols = 1
fig = plt.figure(3)
plt.clf()

ax = plt.subplot(nrows, ncols, 1)
plt.plot(focs, Rlight, label='All optical backscatter')

plt.ylim(0,Rlight.max())
plt.title('All optical backscatter', fontsize=14)
plt.xlabel('Laser focus, distance past column surface ($\mu m$)', fontsize=14)
plt.legend()


ax = plt.subplot(nrows, ncols, 2)
plt.plot(focs, Rhalf, label='Half omega backscatter')
plt.plot(focs, Romega, label='Omega backscatter')

plt.ylim(0,max(Rhalf.max(), Romega.max()))
plt.title('Low frequency backscatter', fontsize=14)
plt.xlabel('Laser focus, distance past column surface ($\mu m$)', fontsize=14)
plt.legend()


ax = plt.subplot(nrows, ncols, 3)
plt.plot(focs, Rblue, label='Two omega backscatter')
plt.plot(focs, Rgreen, label='Three-halves omega backscatter')

plt.ylim(0,max(Rblue.max(), Rgreen.max()))
plt.title('High frequency backscatter', fontsize=14)
plt.xlabel('Laser focus, distance past column surface ($\mu m$)', fontsize=14)
plt.legend()

plt.tight_layout()


#plt.title("", fontsize=20, fontweight='bold')
#plt.vlines(0, -10, 10, linestyle='--', alpha=0.2)
#plt.ylim(2,9)
#
#ax = plt.subplot(nrows, ncols, 2)
#sizes = effics2/np.mean(effics2)*1200
#plt.scatter((focs+6), plas, s=sizes, c = 'r', alpha=0.5, label='Low-energy efficiency')
#plt.ylabel('Plasma scale length ($\mu m$)', fontsize=14)
#plt.xlabel('Laser focus, distance after critical density ($\mu m$)', fontsize=14)
##plt.legend(loc=2)
#plt.title("Low-energy conversion efficiency, $120<KE<500 keV, \pm 40^{\circ}$", fontsize=20, fontweight='bold')
#plt.vlines(0, -10, 10, linestyle='--', alpha=0.2)
#plt.ylim(2,9)
#
#
#ax = plt.subplot(nrows, ncols, 3)
#sizes = Rlight/np.mean(Rlight)*1200
#plt.scatter((focs+6), plas, s=sizes, c = 'g', alpha=0.5, label='All optical reflectivity')
#plt.ylabel('Plasma scale length ($\mu m$)', fontsize=14)
#plt.xlabel('Laser focus, distance after critical density ($\mu m$)', fontsize=14)
##plt.legend(loc=2)
#plt.title("All backscatter", fontsize=20, fontweight='bold')
#plt.vlines(0, -10, 10, linestyle='--', alpha=0.2)
#plt.ylim(2,9)
#
#ax = plt.subplot(nrows, ncols, 4)
#sizes = Rhalf/np.mean(Rhalf)*1200
#plt.scatter((focs+6), plas, s=sizes, c = 'b', alpha=0.5, label='Half-omega backscatter')
#plt.ylabel('Plasma scale length ($\mu m$)', fontsize=14)
#plt.xlabel('Laser focus, distance after critical density ($\mu m$)', fontsize=14)
##plt.legend(loc=2)
#plt.title("Half-omega backscatter", fontsize=20, fontweight='bold')
#plt.vlines(0, -10, 10, linestyle='--', alpha=0.2)
#plt.ylim(2,9)
#
#ax = plt.subplot(nrows, ncols, 5)
#sizes = Romega/np.mean(Romega)*1200
#plt.scatter((focs+6), plas, s=sizes, c = 'b', alpha=0.5, label='Omega backscatter')
#plt.ylabel('Plasma scale length ($\mu m$)', fontsize=14)
#plt.xlabel('Laser focus, distance after critical density ($\mu m$)', fontsize=14)
##plt.legend(loc=2)
#plt.title("Omega backscatter", fontsize=20, fontweight='bold')
#plt.vlines(0, -10, 10, linestyle='--', alpha=0.2)
#plt.ylim(2,9)
#
#ax = plt.subplot(nrows, ncols, 6)
#sizes = Rblue/np.mean(Rblue)*1200
#plt.scatter((focs+6), plas, s=sizes, c = 'b', alpha=0.5, label='Two omega backscatter')
#plt.ylabel('Plasma scale length ($\mu m$)', fontsize=14)
#plt.xlabel('Laser focus, distance after critical density ($\mu m$)', fontsize=14)
##plt.legend(loc=2)
#plt.title("Two omega backscatter", fontsize=20, fontweight='bold')
#plt.vlines(0, -10, 10, linestyle='--', alpha=0.2)
#plt.ylim(2,9)
#
#ax = plt.subplot(nrows, ncols, 7)
#sizes = Rgreen/np.mean(Rgreen)*1200
#plt.scatter((focs+6), plas, s=sizes, c = 'b', alpha=0.5, label='Two omega backscatter')
#plt.ylabel('Plasma scale length ($\mu m$)', fontsize=14)
#plt.xlabel('Laser focus, distance after critical density ($\mu m$)', fontsize=14)
##plt.legend(loc=2)
#plt.title("3/2 omega backscatter", fontsize=20, fontweight='bold')
#plt.vlines(0, -10, 10, linestyle='--', alpha=0.2)
#plt.ylim(2,9)
#plt.text(30, 6, 'Note: 3$\omega$/2 spectral peaks\nnot seen in most cases', fontsize=18)
#
#plt.tight_layout()
#
##ax.set_ylim(0, 2.5)
#
##fig = plt.figure(2)
##sns.jointplot(focs, effics)
#
##ax.plot(focs, Rlight/10, 'k.', label = 'Light reflectivity / 10')
##ax.plot(plas, Rgreen, 'g.', label = 'Green reflectivity')
##ax.plot(effics, Rhalf/5, 'm.', label = 'Half$\omega$ reflectivity')
##ax.plot(plas, effics*8, 'r.', label = 'Efficiency (120/40)')
#
