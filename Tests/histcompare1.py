# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:45:58 2016

@author: Scott
"""
import gzip
import os
import numpy as np
import lspreader2 as rd
import scipy.interpolate as interp

import matplotlib.pyplot as plt
p4dir1 = r'C:\Users\Scott\Documents\temp\mar1test\fully'
p4dir2 = r'C:\Users\Scott\Documents\temp\mar1test\triply'

values, labels = rd.read_hist(p4dir1)
time1 = values[0]*1e6 # Times in fs
energy1 = values[8]
cputime1 = values[7]

values, labels = rd.read_hist(p4dir2)
time2 = values[0]*1e6 # Times in fs
energy2 = values[8]
cputime2 = values[7]

print labels
fig = plt.figure(1)
plt.clf()
plt.plot(time1, energy1, label=os.path.split(p4dir1)[1])
plt.plot(time2, energy2, label=os.path.split(p4dir2)[1])
plt.xlabel('Time (fs)')
plt.ylabel('Net energy (J)')
plt.title('LSP Probe # 8, Net Energy vs. Sim. time')
plt.legend()
#fig.savefig(r'C:\Users\Scott\Documents\LSP\Presentations\Negative energy.png')


fig = plt.figure(2)
plt.clf()
plt.plot(time1, cputime1, label=os.path.split(p4dir1)[1])
plt.plot(time2, cputime2, label=os.path.split(p4dir2)[1])
plt.xlabel('Time (fs)')
plt.ylabel('CPU time (secs)')
plt.title('LSP Probe # 7, CPU time vs. Sim. time')
plt.legend()
#fig.savefig(r'C:\Users\Scott\Documents\LSP\Presentations\Negative energy.png')

cputime2_interp = interp.griddata(time2, cputime2, time1)

fig = plt.figure(3)
plt.clf()
plt.plot(time1, cputime1/cputime2_interp)
plt.xlabel('Time (fs)')
plt.ylabel('CPU time 1 / CPU time 2')
plt.title('LSP Probe # 7, CPU time ratio vs. Sim. time')


print "Simulation 1 took", np.round(np.sum(cputime1)/60/60, 1), "hours"
print "Simulation 2 took", np.round(np.sum(cputime2)/60/60, 1), "hours"
