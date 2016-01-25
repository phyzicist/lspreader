# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 09:55:09 2016

How to convert a frequency plot into a wavelength plot? Shown here.

@author: Scott
"""
import matplotlib.pyplot as plt
import numpy as np

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    

freq = np.linspace(0, 2.5, 120)
Jfreq = gaussian(freq, 1, 0.1) + gaussian(freq, 1.5, 0.1) + gaussian(freq, 2.0, 0.1)

df = np.mean(np.diff(freq))
Jftot = np.sum(Jfreq) * df
print 'Frequency total energy: ' + str(Jftot)

## Make a chop above 0.2 freqs and below 3 freqs
chp = (freq > 0.2) & (freq < 3.0)
freq = freq[chp]
Jfreq = Jfreq[chp]

wlen_fund = 800
wlen = 1/freq * wlen_fund
Jwlen = freq**2 * Jfreq

fig = plt.figure(1)
plt.clf()
ax = plt.subplot(2,1,1)
ax.plot(freq, Jfreq)
ax.set_title('Frequency plot')
ax.set_xlabel('Frequency (in units of the fundamental)')
ax.set_ylabel('Joules/freq.')

ax = plt.subplot(2,1,2)
ax.plot(wlen, Jwlen)
ax.set_title('Wavlength plot')
ax.set_xlabel('Wavelength (in units of the fundamental)')
ax.set_ylabel('Joules/wlen.')
ax.set_xlim(100, 900)

plt.tight_layout()