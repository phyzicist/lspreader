# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:20:31 2016

LSP particle extraction plane analysis.

Scott's consolidated analysis of pext.p4 files (extractions of electrons from four edges.) Very customized.

@author: Scott
"""

# Matplotlib stuff
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from distutils.version import LooseVersion
if LooseVersion(mpl.__version__) < LooseVersion('1.5.0'):    
    # http://stackoverflow.com/questions/11887762/how-to-compare-version-style-strings and 
    print "Matplotlib", mpl.__version__, "might not have colormap 'viridis'. Importing from local colormaps.py."
    import colormaps as cmaps
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    plt.register_cmap(name='inferno', cmap=cmaps.inferno)
    plt.register_cmap(name='magma', cmap=cmaps.magma)
    plt.register_cmap(name='plasma', cmap=cmaps.plasma)

# Other dependencies
import numpy as np
import angular # Radial histogram plot, from Gregory's local file 'angular.py'
import os
import lspreader2 as rd # Hope this works!
import h5stitch2D as hs

def fullAnalyze(p4dir, outdir = '.', shortname = 'Sim'):
    fns = hs.getfnsp4(p4dir, prefix='pext')
    if len(fns) < 1:
        msg = "No pext*.p4(.gz) files found in the directory: " + p4dir
        raise Exception(msg)
    pextarr = loadPexts(fns)
    plotme(pextarr, outdir=outdir, shortname=shortname) # Plots and saves
    
def loadPexts(fns):
    """ Calls lspreader2 to load all the pext*.p4 from a given list of filenames, and concatenates into a single array
    Inputs:
        fns: List of filenames
    Outputs:
        pextarr: A NumPy array of the same form as the output of lspreader2.read_pext(), containing a concatenation of all the pext file contents as opposed to just one
    """
    pextarr = None # An array which will hold the contents of all the pext files
    for i in range(len(fns)):
        fn = fns[i] # Current filename
        pextarr_tmp = rd.read_pext2(fn)
        if i == 0:
            pextarr = pextarr_tmp # Initialize the output array 'outs'
        else:
            pextarr = np.concatenate((pextarr, pextarr_tmp), 0) # Append this pext*.p4 'out' to the output array 'outs'
    return pextarr

def plotme(pextarr, outdir='.', shortname='Sim', mksub=False):
    """
    Make Scott's preferred set of plots (and save to file) from the output (or concatenated outputs) of lspreader "pext()" of pext.
    Inputs:
        pextarr: NumPy array, of the fancy form output by pextarr = lspreader2.read_pext('pext1.p4'). Several can be concatenated as desired.
        outdir: string, path to a folder (already existent) that can hold the outputs
        shortname: string, a very short identifier for the simulation which will go in all the titles and filenames
        mksub: bool, if True then make a subdirectory in 'outdir', otherwise just dump the PNGs directly into 'outdir'
    Outputs:
        Will save PNG figures summarizing the extraction planes into the folder specified by 'outdir'.
    
    TODO: stretch out the time of flight figure
    TODO: fix up the time of flight figure so it's not such an ugly, long string of plotting functions (copied and pasted text)
    """

    # Make a subdirectory for the PNG files, if mksub was set to True
    if mksub:
        plotdir = subdir(outdir, 'Extraction planes')
    else:
        plotdir = outdir
        
    ## GET PHI FROM LSPREADER OUTPUTS OF PEXT (Concatenated into one, of course.)
    massE = 0.511 # electron rest mass, in MeV
    d = {}
    d['t'] = pextarr['t']*1e6 # Time, in fs
    u_norm = np.sqrt(pextarr['ux']**2+pextarr['uy']**2+pextarr['uz']**2) # 'uz', etc. are actually momenta; in what they call gamma-beta units, so unitless and related only to speed of light
    d['KE'] =(np.sqrt(u_norm**2+1)-1)*massE # Electron kinetic energy, in MeV
    d['q'] = -pextarr['q']*1e3 # Amount of negative charge, nanoColoumbs/cm
    d['phi'] = np.arctan2(pextarr['uz'], pextarr['ux'])
        
    ## RADIAL PLOT WITHOUT CUTOFFS
    fig = plt.figure(3)
    plt.clf()
    title=r'Radial electron spectrum'
    angular.angular2(d['q'], d['phi'], d['KE'], fig = fig, cmap='viridis', minQ=0, maxQ=1.5, title=title)
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Radial.png'))
   
    ## RADIAL PLOT WITH CUTOFF
    fig = plt.figure(2)
    condit2 =np.logical_and(d['KE'] > .300, np.abs(d['phi']) > np.deg2rad(180 - 40))
    plt.clf()
    #title = shortname + r": 300/40 charge: " + str(round(np.sum(d['q'][condit2])))
    title=r'Radial electron spectrum, $\pm$40$^{\circ}$'
    angular.angular2(d['q'][condit2], d['phi'][condit2], d['KE'][condit2], fig = fig, cmap='viridis', minQ=0, maxQ=1.5, title=title)
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Radial with cutoffs.png'))
    
    ## BACKWARD ELECTRON SPECTRUM LINE PLOT
    fig = plt.figure(4)
    plt.clf()
    cdt_inner = np.abs(d['phi']) > np.deg2rad(180 - 6.3)
    cdt_outer = np.abs(d['phi']) > np.deg2rad(180 - 40)
    nbins = 60
    Erange = np.array([0, 3.0])
    dE = (np.max(Erange) - np.min(Erange))/nbins # Width, in MeV, of each histogram bin
    [histvals, histedges] = np.histogram(d['KE'][cdt_outer], range=Erange, bins=nbins, weights=d['q'][cdt_outer])
    picY1 = histvals / dE # convert from charge/bin to nC/MeV/cm
    picX1 = histedges[1:] - (histedges[1] - histedges[0])/2
    [histvals, histedges] = np.histogram(d['KE'][cdt_inner], range=Erange, bins=nbins, weights=d['q'][cdt_inner])
    picY2 = histvals / dE # convert from charge/bin to nC/MeV/cm
    picX2 = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picX1, picY1, picX2, picY2)
    plt.xlabel("Kinetic energy, MeV")
    plt.ylabel("Charge density, 2D (nC/MeV/cm)")
    plt.ylim(0, 625)
    plt.legend([r'Backward $\pm$40$^{\circ}$', r'Backward $\pm$6.3$^{\circ}$'])
    plt.title(r'Electron spectrum', fontsize=20)
    plt.tight_layout()
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Electron spectrum.png'))

    
    ## TIME OF FLIGHT PLOT
    fig = plt.figure(5)
    tmax = np.max(d['t']) # The time of the end of the sim, presumably. Draw a vertical line at this time.
    trange = (100,400) # The plotting range
    plt.clf()
    
    ax = plt.subplot(411)
    cdt_time = np.logical_and(cdt_outer, d['KE'] > 3.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    plt.title('KE > 3 MeV')
    
    
    ax = plt.subplot(412)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 3.0, d['KE'] > 1.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    plt.title('1 MeV < KE < 3 MeV')
    
    ax=  plt.subplot(413)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 1.0, d['KE'] > .5)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    plt.title('500 keV < KE < 1 MeV')
    
    ax = plt.subplot(414)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < .5, d['KE'] > .1)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    plt.xlabel("Time (fs)")
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    plt.title('100 keV < KE < 500 keV')
    
    plt.suptitle('Electron time of flight, $\pm$40$^{\circ}$', fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Electron time of flight.png'))


    ## TIME OF FLIGHT PLOT, Zoomed in
    fig = plt.figure(6)
    tmax = np.max(d['t']) # The time of the end of the sim, presumably. Draw a vertical line at this time.
    tmid = np.sum(d['t']*d['q'])/np.sum(d['q']) # Centroid of charge deposition, in time
    trange = (tmid - 50,tmid + 50) # The plotting range
    plt.clf()
    
    ax = plt.subplot(411)
    cdt_time = np.logical_and(cdt_outer, d['KE'] > 3.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('KE > 3 MeV')
    
    
    ax = plt.subplot(412)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 3.0, d['KE'] > 1.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('1 MeV < KE < 3 MeV')
    
    ax=  plt.subplot(413)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 1.0, d['KE'] > .5)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('500 keV < KE < 1 MeV')
    
    ax = plt.subplot(414)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < .5, d['KE'] > .1)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    plt.xlabel("Time (fs)")
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('100 keV < KE < 500 keV')
    
    plt.suptitle('Electron time of flight (ZOOMED), $\pm$40$^{\circ}$', fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Zoomed TOF.png'))


def subdir(folder, name):
    """ Make a subdirectory in the specified folder, if it doesn't already exist"""
    subpath = os.path.join(folder,name)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    return subpath