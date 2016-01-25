# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 11:16:47 2016

Secondary analysis.

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
plas = []
for name in namelist:
    m1 = re.search('^a([0-9]*?)fp{0,1}(-*[0-9]*?)_mres_so', name)
    if m1:
        plas.append(float(m1.groups()[0])*0.1)
        focs.append(float(m1.groups()[1]))
        goodnames.append(name)
    
print goodnames
print focs
print plas

effics = []
Rlight = []
Rgreen = []
Rhalf = []

for name in goodnames:
    fn_effic = os.path.join(analydir, name, name + ' - Electron Efficiency.csv')
    arr_effic = np.genfromtxt(fn_effic, skip_header=2, delimiter=',')
    effics.append(arr_effic[0,1])
    
    fn_refl = os.path.join(analydir, name, name + ' - Reflectivity.csv')
    arr_refl = np.genfromtxt(fn_refl, skip_header=1, delimiter=',')
    Rlight.append(arr_refl[0,1] - arr_refl[2,1]) # All reflectivity energy, less the static
    Rgreen.append(arr_refl[3,1])
    Rhalf.append(arr_refl[6,1])


focs = np.array(focs)
effics = np.array(effics)
Rlight = np.array(Rlight)
Rgreen = np.array(Rgreen)
Rhalf = np.array(Rhalf)
plas = np.array(plas)

print effics
print Rlight

fig = plt.figure(1)
plt.clf()
ax = plt.subplot(111)
#ax.plot(focs, Rlight/10, 'k.', label = 'Light reflectivity / 10')
#ax.plot(plas, Rgreen, 'g.', label = 'Green reflectivity')
ax.plot(effics, Rhalf/5, 'm.', label = 'Half$\omega$ reflectivity')
#ax.plot(plas, effics*8, 'r.', label = 'Efficiency (120/40)')
plt.legend()
#ax.set_ylim(0, 2.5)
