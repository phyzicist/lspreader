# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 09:55:09 2016

How to convert a frequency plot into a wavelength plot? Shown here.

Also, adding in the effect of a TLMB mirror by loading its transmission from CSV.

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

import matplotlib.pyplot as plt
import numpy as np
from special import tlmb

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

freq = np.linspace(0, 5, 400)
Jfreq = 10*gaussian(freq, 1, 0.1) + gaussian(freq, 1.5, 0.1) + gaussian(freq, 2.0, 0.1) + gaussian(freq, 0.5, 0.05)

df = np.mean(np.diff(freq))
Jftot = np.sum(Jfreq) * df
print 'Frequency total energy: ' + str(Jftot)

## Make a chop above 0.1 (because below that is essentially electrostatic)
chp = (freq > 0.1) # Chopping condition
wlen_fund = 800
wlen = (1/freq[chp] * wlen_fund)[::-1] # The -1 is to reverse the arrays
Jwlen = (freq[chp]**2 * Jfreq[chp] / wlen_fund)[::-1]

## Adjust by the TLMB transmission
T_tlmb = tlmb.trans(wlen)
Jwlen2 = Jwlen * T_tlmb

fig = plt.figure(1)
plt.clf()
ax = plt.subplot(2,1,1)
ax.plot(freq, Jfreq)
ax.set_title('Frequency plot')
ax.set_xlabel('Frequency (in units of the fundamental)')
ax.set_ylabel('Joules/freq.')

ax = plt.subplot(2,1,2)
ax2 = ax.twinx()
ax2.plot(wlen, T_tlmb, 'g')
ax2.set_ylim(0, 1)
ax.plot(wlen, Jwlen, 'b')
ax.plot(wlen, Jwlen2, 'k')
ax.set_title('Wavlength plot')
ax.set_xlabel('Wavelength (in units of the fundamental)')
ax.set_ylabel('Joules/wlen.')
ax.set_xlim(100, 1200)

fig = plt.figure(2)
plt.clf()
ax = plt.subplot(2,1,1)
ax.set_title('Frequency plot')
ax.set_xlabel('Frequency (in units of the fundamental)')
ax.set_ylabel('Joules/freq.')
ax.plot(wlen, Jwlen2, 'k')
ax.set_xlim(100, 900)

fig = plt.figure(2)
plt.clf()
ymin = 0
ymax = np.max(Jwlen2[(wlen > 100) & (wlen < 600)])
ax = plt.subplot(2,1,1)
ax.set_title('Blue-green light')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Joules/ optical nm')
ax.plot(wlen, Jwlen2, 'r')
ax.set_xlim(200, 900)
ax.vlines([wlen_fund*2, wlen_fund, wlen_fund * 2./3., wlen_fund / 2.0], ymin, ymax, colors=[(0.6,0,0),'r','g','b'], linestyle=':') # vertical lines, colored by frequency
ax.set_ylim(0, ymax)

ax = plt.subplot(2,1,2)
ymin = 0
ymax = np.max(Jwlen2[(wlen > 950) & (wlen < 1700)])
ax.set_title('Mid-IR light')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Joules/ optical nm')
ax.plot(wlen, Jwlen2, 'r')
ax.set_xlim(900, 1800)
ax.vlines([wlen_fund*2, wlen_fund, wlen_fund * 2./3., wlen_fund / 2.0], ymin, ymax, colors=[(0.6,0,0),'r','g','b'], linestyle=':') # vertical lines, colored by frequency
ax.set_ylim(ymin, ymax)
plt.tight_layout()