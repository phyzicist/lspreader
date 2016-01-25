# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 09:55:09 2016

How to convert a frequency plot into a wavelength plot? Shown here.

Also, adding in the effect of a TLMB mirror by loading its transmission from CSV.

And now, testing out the possibility of running it in scottplots.

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

def wlenPlot(freq, Jfreq, wlen_fund=800):
    """ Make the backscatter wavelength plot from data for a frequency plot. Don't save to file, just return the handle.
    Inputs:
        freq: 1D array, list of frequencies (in units of the fundamental) (the X axis for a power spectrum plot)
        Jfreq: 1D array, same length as freq, list of Joules/freq. at those frequencies (the Y axis for a power spectrum plot)
        wlen_fund: (optional) number, wavelength of the fundamental, in nm
    Outputs:
        fig: handle to the plot where X axis is wavelength (in nm), and Y axis is Joules/wlen. (which can be saved to file, e.g. with fig.savefig())
    """
    ## Make a chop above 0.1 (because below that is essentially electrostatic), and then transform the array
    chp = (freq > 0.1) # Chopping condition
    wlen = (1/freq[chp] * wlen_fund)[::-1] # wlen = 1 / frequency. The -1 is to reverse the arrays, since wavelength and frequency are reversed
    Jwlen = (freq[chp]**2 * Jfreq[chp] / wlen_fund)[::-1] # Adjust the Y axis, based on the fact that dwlen/dfreq = -1/freq^2
    
    ## Adjust by the TLMB transmission
    T_tlmb = tlmb.trans(wlen)
    Jwlen2 = Jwlen * T_tlmb
    
    fig = plt.figure(2)
    plt.clf()
    ymin = 0
    ymax = np.max(Jwlen2[(wlen > 100) & (wlen < 600)])
    ax = plt.subplot(2,1,1)
    ax.set_title('Blue-green backscatter through TLMB')
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Joules/ optical nm')
    ax.plot(wlen, Jwlen2, 'b')
    ax.set_xlim(200, 900)
    ax.vlines([wlen_fund*2, wlen_fund, wlen_fund * 2./3., wlen_fund / 2.0], ymin, ymax, colors=[(0.6,0,0),'r','g','b'], linestyle=':') # vertical lines, colored by frequency
    ax.set_ylim(0, ymax)
    
    ax = plt.subplot(2,1,2)
    ymin = 0
    ymax = np.max(Jwlen2[(wlen > 950) & (wlen < 1700)])
    ax.set_title('Mid-IR backscatter through TLMB')
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Joules/ optical nm')
    ax.plot(wlen, Jwlen2, 'b')
    ax.set_xlim(900, 1800)
    ax.vlines([wlen_fund*2, wlen_fund, wlen_fund * 2./3., wlen_fund / 2.0], ymin, ymax, colors=[(0.6,0,0),'r','g','b'], linestyle=':') # vertical lines, colored by frequency
    ax.set_ylim(ymin, ymax)
    plt.tight_layout()
    
    return fig
    

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

freq = np.linspace(0, 5, 2000)
Jfreq = 10*gaussian(freq, 1, 0.1) + gaussian(freq, 1.5, 0.1) + gaussian(freq, 2.0, 0.1) + gaussian(freq, 0.5, 0.05)

df = np.mean(np.diff(freq))
Jftot = np.sum(Jfreq) * df
print 'Frequency total energy: ' + str(Jftot)
fig = wlenPlot(freq, Jfreq)