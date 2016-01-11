# -*- coding: utf-8 -*-
"""
Analyze a batch of p4 files for frequency components

Created on Wed Dec 30 15:09:37 2015

@author: Scott
"""
import h5py
import os
import lspreader2 as rd
import numpy as np
from h5stitch2D import chunkIt, getfnsp4, fields2D, h5fields2D, h5fields2Dpar2, getTimes
import gzip
from copy import deepcopy

# Matplotlib stuff
import matplotlib as mpl
mpl.use('Agg')
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

try:
    from mpi4py import MPI
except:
    print "WARNING: MPI4PY FAILED TO LOAD. DO NOT CALL PARALLEL FUNCTIONS."
    
    
def freqanalyze(fns, data=None, datpath=None, fld_id = 'Ez', divsp=2):
    # Inputs:
    #   divsp: integer, divisor by which to reduce the spatial resolution (e.g. divsp = 2 reduces field dimensions from 300x200 to 150x100)
    #   data: (Optional) A python dictionary returned by h5stitch2D.py "fields2D()". If none is supplied, it will be read in by reading from datah5 (fields2D.hdf5) or by looking at p4 filenames.
    #   datpath: (Optional) A string path of the hdf5 file from which data can be loaded.
    #   fns: list of filenames to analyze, if neither data nor datah5 are given.
    #   h5path 
    #   fld_id: A string specifying the field component to analyze. This field must be contained within the "data" dictionary
    
    if data:
        print "Data supplied at input to freqanalyze(), so no need to read in the files. Moving on to frequency analysis."
        times = data['times']*1e6 # times, converted to fs
        fns = data['filenames']
        xgv = data['xgv'][::divsp]*1e4 # spatial X, converted to microns
        zgv = data['zgv'][::divsp]*1e4 # spatial Z, converted to microns
        Ez = data[fld_id][:,::divsp,::divsp] ## I call it "Ez" as a variable name, but this could be any field.
    elif datpath:
        print "Will read fields2D data from the HDF5 file."
        with h5py.File(datpath, 'r') as data:
            times = data['times'][...]*1e6 # times, converted to fs
            fns = data['filenames'][...]
            xgv = data['xgv'][::divsp][...]*1e4 # spatial X, converted to microns
            zgv = data['zgv'][::divsp][...]*1e4 # spatial Z, converted to microns
            Ez = data[fld_id][:,::divsp,::divsp][...] ## I call it "Ez" as a variable name, but this could be any field.
    else:
        print "Will read ", len(fns), "files."
        # Load in the data
        data = fields2D(fns, fld_ids = [fld_id])
        #print "Data read in."
        times = data['times']*1e6 # times, converted to fs
        fns = data['filenames']
        xgv = data['xgv'][::divsp]*1e4 # spatial X, converted to microns
        zgv = data['zgv'][::divsp]*1e4 # spatial Z, converted to microns
        Ez = data[fld_id][:,::divsp,::divsp] ## I call it "Ez" as a variable name, but this could be any field.

    data2 = {}
    data2['times_fs'] = times
    data2['xgv_um'] = xgv
    data2['zgv_um'] = zgv
    data2['fns'] = fns

    #print "times, fns, xgv, zgv, Ez"
    #print times.shape
    #print fns.shape
    #print xgv.shape
    #print zgv.shape
    #print Ez.shape

    # Assume equal spacing in time and space, and get the deltas
    #print "Calculating dt, dx, dz"
    dt = np.mean(np.diff(times))
    dx = np.mean(np.diff(xgv))
    dz = np.mean(np.diff(zgv))

    #print dt, dx, dz

    # Calculate the frequency of the laser from its wavelength (800 nm)
    c = 3e8 # Speed of light in m/s
    wl = 0.8e-6 # Wavelength of laser in m
    fr_fund = c/wl # Frequency of the laser (the 'fundamental'), in Hz

    freqHz = np.fft.rfftfreq(len(times), d = dt/1e15) # Frequencies in Hz, for upcoming real FFT
    freq = freqHz/fr_fund # Frequencies in units of the fundamental

    data2['freq'] = freq

    # Make some interesting maps
    Imap = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map
    FTmap1 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
    FTmap2 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
    FTmap3 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
    FTmap4 = np.zeros((Ez.shape[1],Ez.shape[2])) # An intensity map, but only over certain frequencies
    pwr_sum = np.zeros(freq.shape)

    map1_condit = np.logical_and(freq > 0.3, freq < 0.7)
    map2_condit = np.logical_and(freq > 0.9, freq < 1.1)
    map3_condit = np.logical_and(freq > 1.3, freq < 1.7)
    map4_condit = np.logical_and(freq > 1.8, freq < 2.3)

    #print "Doing calculations"
    for i in range(len(zgv)):
        for j in range(len(xgv)):
            tline = Ez[:,i,j]
            
            # Make a simple intensity map        
            esum = np.sum(tline**2)
            Imap[i,j] = esum
            
            # Make a more complex ftransform map
            time_ft = np.fft.rfft(tline)
            pwr = np.abs(time_ft)**2
            FTmap1[i,j] = np.sum(pwr[map1_condit])
            FTmap2[i,j] = np.sum(pwr[map2_condit])
            FTmap3[i,j] = np.sum(pwr[map3_condit])
            FTmap4[i,j] = np.sum(pwr[map4_condit])
            pwr_sum = pwr_sum + pwr;



    data2['pwr_sum'] = pwr_sum
    data2['Imap'] = Imap
    data2['FTmap_0_5'] = FTmap1
    data2['FTmap_1'] = FTmap2
    data2['FTmap_1_5'] = FTmap3
    data2['FTmap_2'] = FTmap4

    return data2

def plotme(data2, folder='', pltdict = None):
    xgv = data2['xgv_um']
    zgv = data2['zgv_um']
    freq = data2['freq']
    times = data2['times_fs']
    pwr_sum = data2['pwr_sum']
    Imap = data2['Imap']
    FTmap1 = data2['FTmap_0_5']
    FTmap2 = data2['FTmap_1']
    FTmap3 = data2['FTmap_1_5']
    FTmap4 = data2['FTmap_2']

    meantime = np.mean(times)
    tplus = np.max(times) - np.mean(times)
    tstring = 't=' + "{:.1f}".format(meantime) + " fs $\pm$ " + "{:.1f}".format(tplus) + " fs"
    
    # Make some plots
    #print "Making plots."
    figs = []

    # Time-integrated intensity map (a.u.)
    fig = plt.figure(1)
    figs.append(fig)
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    im = ax.pcolorfast(xgv,zgv,Imap, cmap='gray')
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Z (um)')
    ax.set_title('Integrated power from Ez (a.u.)')
    X,Z = np.meshgrid(xgv,zgv)
    plt.contour(X,Z,Imap, cmap='viridis')
    ax.text(-19, 17, tstring, fontsize=24, color='white')
    if pltdict:
        cmin = pltdict['Imap']['min']
        cmax = pltdict['Imap']['max']
        im.set_clim(vmin=cmin, vmax=cmax)

    # Power spectrum
    fig = plt.figure(2)
    figs.append(fig) # Add the figure to the figures list
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    ax.plot(freq, np.log10(pwr_sum))
    ax.set_xlabel('Frequency (normalized to laser fundamental)')
    ax.set_ylabel('Log10(Power spectrum) of Ez (a.u.)')
    ax.set_title('Log of Power spectrum at ' + tstring)
    ax.set_xlim(0, 2.5)
    ax.xaxis.grid() # vertical lines
    if pltdict:
        ymin = np.log10(pltdict['pwr_sum']['min'])
        ymax = np.log10(pltdict['pwr_sum']['max'])
        ax.set_ylim(ymin, ymax)
    
    # Fourier transform subset map 1
    fig = plt.figure(3)
    figs.append(fig)
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    im = ax.pcolorfast(xgv,zgv,FTmap1, cmap='viridis')
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Z (um)')
    ax.set_title('Ez, Half omega map (0.3 to 0.7 x fundamental) (a.u.)')
    ax.text(-19, 17, tstring, fontsize=24, color='white')
    fig.colorbar(im, label='Frequency-cut power spectrum integral (a.u.)')
    if pltdict:
        cmin = pltdict['FTmap_0_5']['min']
        cmax = pltdict['FTmap_0_5']['max']
        im.set_clim(vmin=cmin, vmax=cmax)

    # Fourier transform subset map 2
    fig = plt.figure(4)
    figs.append(fig)
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    im = ax.pcolorfast(xgv,zgv,FTmap2, cmap='viridis')
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Z (um)')
    ax.set_title('Ez, Fundamental map (0.9 to 1.1 x fundamental) (a.u.)')
    ax.text(-19, 17, tstring, fontsize=24, color='white')
    fig.colorbar(im, label='Frequency-cut power spectrum integral (a.u.)')
    if pltdict:
        cmin = pltdict['FTmap_1']['min']
        cmax = pltdict['FTmap_1']['max']
        im.set_clim(vmin=cmin, vmax=cmax)

    # Fourier transform subset map 3
    fig = plt.figure(5)
    figs.append(fig)
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    im = ax.pcolorfast(xgv,zgv,FTmap3, cmap='viridis')
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Z (um)')
    ax.set_title('Ez, Three-halves omega map (1.3 to 1.7 x fundamental) (a.u.)')
    ax.text(-19, 17, tstring, fontsize=24, color='white')
    fig.colorbar(im, label='Frequency-cut power spectrum integral (a.u.)')
    if pltdict:
        cmin = pltdict['FTmap_1_5']['min']
        cmax = pltdict['FTmap_1_5']['max']
        im.set_clim(vmin=cmin, vmax=cmax)

    # Fourier transform subset map 4
    fig = plt.figure(6)
    figs.append(fig)
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    im = ax.pcolorfast(xgv,zgv,FTmap4, cmap='viridis')
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Z (um)')
    ax.set_title('Ez, Two omega map (1.8 to 2.3 x fundamental) (a.u.)')
    ax.text(-19, 17, tstring, fontsize=24, color='white')
    fig.colorbar(im, label='Frequency-cut power spectrum integral (a.u.)')
    if pltdict:
        cmin = pltdict['FTmap_2']['min']
        cmax = pltdict['FTmap_2']['max']
        im.set_clim(vmin=cmin, vmax=cmax)

    print "Saving figures"
    for i in range(len(figs)):
        fn = os.path.join(folder, 'fig' + str(i) + '.png')
        figs[i].savefig(fn)

    return figs

def saveData2(data2, folder = '', name = 'data2.hdf5'):
    """Save the data2 dictionary to an HDF5 file, assuming it contains only NumPy arrays"""
    print "Saving HDF5 file."
    h5path = os.path.join(folder, name)
    with h5py.File(h5path, 'w') as f:
        for k in data2:
            f.create_dataset(k, data = data2[k], compression='gzip', compression_opts=4)

    return h5path

def loadData2(h5path):
    """Read the data2 HDF5 file back into the data2 array"""
    print "Reading in HDF5 file."
    data2 = {}
    with h5py.File(h5path, 'r') as f:
        for k in f.keys():
            data2[k] = f[k][...]
    return data2

def chunkFx(mylist, n):
    """Break list into fixed, n-sized chunks. The final element of the new list will be n-sized or less"""
    # Modified from http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    chunks = []
    for i in range(0, len(mylist), n):
        chunks.append(mylist[i:i+n])
    return chunks

def unChunk(l):
    """Flatten the first dimension of a list. E.g. if input is l = [[1,2,],[3,4]], output is [1,2,3,4]"""
    # Copied from http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    return [item for sublist in l for item in sublist]

def getPltDict(data2):
    """ Get the plotting dictionary, which includes things such as maximums needed to make plots, from a data2 dictionary"""
    pltdict = {}
    for k in data2:
        pltdict[k] = {}
        # Get the maximum values
        try:
            pltdict[k]['max'] = np.max(data2[k])
        except: # Perhaps it is non-numeric? Bad coding, here.
            pltdict[k]['max'] = np.float('nan')
        
        # Get the minimum values
        try:
            pltdict[k]['min'] = np.min(data2[k])
        except: # Perhaps it is non-numeric? Bad coding, here.
            pltdict[k]['min'] = np.float('nan')

    return pltdict

def bestPltDict(pltdicts):
    """Create the best plot dictionary from a list of plot dictionaries"""
    pltdict = pltdicts[0] # Define starting point for master plot dict. Note that this is not a deep copy, but that's ok in this case.
    for k in pltdict:
        for i in range(len(pltdicts)):
            # Get the widest max and min values possible (compare current plot dict to master plot dict)
            pltdict[k]['max'] = max(pltdicts[i][k]['max'], pltdict[k]['max'])
            pltdict[k]['min'] = min(pltdicts[i][k]['min'], pltdict[k]['min'])
    return pltdict

def mainPar():
    p4dir = r'/tmp/ngirmang.1-2d-nosolid-151113/' # This folder contains the p4 files
    outdir = r'/home/feister.7/lsp/runs/greg_run' # This folder already exists and will store the files in subfolders
    
    nbatch = 45 # Number of files per batch (Each batch is analyzed in Fourier space)
    
    nprocs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()
    comm = MPI.COMM_WORLD
    
    # Get files in folder
    fns = getfnsp4(p4dir)[::2]
    # Split the filenames list into batches of size nbatch (each batch will be fourier-analyzed)
    fns_batched = chunkFx(fns, nbatch)
    # Make some sub-batches
    alt_batched = chunkFx(fns[int(round(nbatch/2)):], nbatch) # Same as above, but offset
    fns_batched = fns_batched + alt_batched
    
    if rank == 0:
        print "NUMBER OF CHUNKS: ", len(fns_batched)

    # Spread the batches of filenames across the various processors
    batches_chunked = chunkIt(fns_batched,nprocs)
    batches_part = batches_chunked[rank] # This processor's chunk of filename batches
    
    print "Processor", rank, "assigned", len(batches_part), "batches."
    
    pltdicts_part = []
    data2s_part = []
    outsubdirs_part = []
    for fns_batch in batches_part:
        data2 = freqanalyze(fns_batch, divsp=1)
        
        # Give the subdirectory for this dataset a name, and make the directory.
        meantime = np.mean(data2['times_fs']) # Mean time, in fs
        dirname = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the folder label be '00512' for t = 51.2342 fs
        outsubdir = os.path.join(outdir, dirname)
        outsubdirs_part.append(outsubdir)

        if not os.path.exists(outsubdir):
            os.mkdir(outsubdir)
        
        # Get plot helper dictionary, which includes maximum values, and will be needed for plotting
        pltdict = getPltDict(data2)
        pltdicts_part.append(pltdict)

        # Save the figures and HDF5 into the subdirectory
        h5path = saveData2(data2, folder=outsubdir)
        
        # Append to the data2 list
        data2s_part.append(data2)

    # Compare all the plot dictionaries, and select the actual plotting parameters
    pltdicts_chunked = comm.gather(pltdicts_part, root=0)
    if rank == 0:
        pltdicts_all = unChunk(pltdicts_chunked)
        print pltdicts_all
        pltdict_final = bestPltDict(pltdicts_all)
        print "FINAL PLOT DICTIONARY: ", pltdict_final
    else:
        pltdicts_all = None
        pltdict_final = None

    pltdict_final = comm.bcast(pltdict_final, root=0)

    # Re-iterate over the chunks for this processor, saving images.

    for i in range(len(outsubdirs_part)):
        data2 = data2s_part[i]
        outsubdir = outsubdirs_part[i]
        plotme(data2, folder=outsubdir, pltdict = pltdict_final)

def totalSer():
    p4dir = r'/tmp/ngirmang.1-2d-nosolid-151113/' # This folder contains the p4 files
    outdir = r'/home/feister.7/lsp/runs/greg_run' # This folder already exists and will store the files in subfolders
    fns = getfnsp4(p4dir)
    data2 = freqanalyze(fns, divsp=1)

    # Give the subdirectory for this dataset a name, and make the directory.
    meantime = np.mean(data2['times_fs']) # Mean time, in fs
    dirname = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the folder label be '00512' for t = 51.2342 fs
    outsubdir = os.path.join(outdir, dirname)
    if not os.path.exists(outsubdir):
        os.mkdir(outsubdir)

    pltdict = getPltDict(data2)
    h5path = saveData2(data2, folder=outsubdir)
    plotme(data2, folder=outsubdir, pltdict=pltdict)

if __name__ == "__main__":
    ## MAIN PROGRAM
    #mainPar() # Frequency analysis of the folder, in multiple steps
    totalSer() # Freqency analysis of all data as a whole

