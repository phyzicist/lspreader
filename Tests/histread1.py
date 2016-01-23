# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:45:58 2016

@author: Scott
"""
import gzip
import os
import numpy as np
import lspreader2 as rd

import matplotlib.pyplot as plt
p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test2'

values, labels = rd.read_hist(p4dir)
time = values[0]*1e6 # Times in fs
energy = values[8]

print labels
fig = plt.figure(1)
plt.clf()
plt.plot(time, energy)
plt.xlabel('Time (fs)')
plt.ylabel('Net energy (J)')
plt.title('LSP Probe # 8, Net Energy vs. Time')
#fig.savefig(r'C:\Users\Scott\Documents\LSP\Presentations\Negative energy.png')