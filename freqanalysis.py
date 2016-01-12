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
import timeit
import shutil

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


# Parallel processing stuff
from multiprocessing import Pool
try:
    from mpi4py import MPI
except:
    print "WARNING: MPI4PY FAILED TO LOAD. DO NOT CALL MPI-BASED FUNCTIONS."
    
# Done with importing modules! Let the games begin.

def freqFull(p4dir, outdir = '', nbatch = 20, divsp = 1, fld_ids = ['Ez', 'Ex', 'By'], npool = 1):
    """ Perform frequency analysis of Ez, Ex, and By (and/or beyond) on many batches of .p4 or .p4.gz; that is, on an entire folder. Analysis can be in parallel or serial."""
    # Inputs:
    #   p4dir: folder containing .p4 and/or .p4.gz field files
    #   outdir: the full path of where the .png plots and their .hdf5 data will be saved
    #   divsp: an integer factor by which to decimate the spatial resolution of the field arrays.
    #   fld_ids: a list of field identifiers (as seen by readflds()) upon which fields frequency analysis will be done
    #   npool: integer number of threads desired in the parallel pool. If npool <= 1, analysis is done in serial
    #   nbatch: the number of files upon which each frequency analysis is performed. Smaller gives better time resolution, but larger gives better frequency resolution.
    # Outputs:
    #   .png and .hdf5 files (with which revised plots can be made) are saved in subdirectories of "outdir"

    print "BEGINNING FREQUENCY ANALYSIS OF DATASET."
    # Start a parallel pool of processors. If npool is 1, run the whole analysis in serial. 
    if npool > 1:
        print "Starting parallel pool of", npool, "threads for this analysis."
        pool = Pool(npool)
    else:
        print "Serial analysis requested (so no pool started)."
        pool = None

    # Look at the files in the folder, and split them into good-sized batches for analysis
    fns_all = getfnsp4(p4dir)

    # Split the filenames list into batches of size nbatch (each batch will be fourier-analyzed)
    main_batched = chunkFx(fns_all, nbatch)
    # Make some sub-batches
    alt_batched = chunkFx(fns_all[int(round(nbatch/2)):], nbatch) # Same as above, but offset
    fns_batched = main_batched + alt_batched


    print "Split files into", len(fns_batched), "small batches of", nbatch, "files and one large batch of", len(fns_all)
    
    ## Analyze the ALL TIMES batch, making plots
    print "Starting with the 'All times' frequency analysis."
    freqBatch(fns_all,outdir = outdir, divsp = divsp,  fld_ids = fld_ids, pool = pool, alltime=True) # Will make the plots as well in this case

    ## Analyze all the TIME-RESOLVED batches, and make plots with consistent colorbar limits
    # Pre-allocate data structures  
    pltdicts = {} # A dictionary with keys "fld_ids", and each key unlocks a list of plot dictionaries (one for each file)
    for fld_id in fld_ids:
        pltdicts[fld_id] = [None]*len(fns_batched) # Pre-allocate the lists of plot dictionaries

    data2s_dict = [None]*len(fns_batched) # A list with each element being a dictionary, and each dictionary having keys "fld_ids", and each key unlocking a data2 dict (which holds the entire result of frequency analysis). Could this become a memory problem? Perhaps.

    for i in range(len(fns_batched)):
        print "Small batch frequency analysis", i, "of", len(fns_batched)
        data2s_dict[i], pltdict_dict = freqBatch(fns_batched[i], outdir = outdir, divsp = divsp,  fld_ids = fld_ids, pool = pool)
        for fld_id in pltdict_dict:
            pltdicts[fld_id][i] = pltdict_dict[fld_id] # Fill in the i-th element in this field's list of plot dictionaries

    # Extract normalization max and min between all time-resolved frequency analyses, for consistent plots across all files
    pltdict = {} # A dictionary with keys "fld_ids", and each key unlocks a single plot dictionary (the best max/min of all files)
    for fld_id in fld_ids:
        pltdict[fld_id] = bestPltDict(pltdicts[fld_id]) # Get the best max and min for plotting

    # Make plots for the time-resolved frequency analyses and save them to file
    for i in range(len(fns_batched)):
        for fld_id in fld_ids:    
            plotme(data2s_dict[i][fld_id], outdir=outdir, pltdict=pltdict[fld_id], fld_id = fld_id) # Make plots and save to png
    plt.close('all') # We are done with plots, so close any that are still open."

    return outdir


def freqBatch(fns, outdir = '', divsp = 1, fld_ids = ['Ez', 'Ex', 'By'], pool = None, alltime=False):
    """ Perform frequency analysis of Ez, Ex, and By (and/or beyond) on a single batch of .p4 or .p4.gz files. Does not make plots, unless alltime = True. Does save the data to .hdf5 files."""
    data2_dict = {} # A dictionary with fld_ids as keys; and each key unlocks its own respective data2 dictionary.
    pltdict_dict = {} # A dictionary with fld_ids as keys; and each key unlocks its own respective pltdict dictionary.
 
    data = fields2D(fns, fld_ids = fld_ids, pool = pool)
    for fld_id in fld_ids:
        data2 = freqanalyze(fns, fld_id = fld_id, data = data, pool = pool, divsp = divsp)
        pltdict = getPltDict(data2) # Get the maxima/minima across multiple files, for homogeneous plotting across time steps
        h5path = freqSave(data2, outdir = outdir, fld_id = fld_id, alltime=alltime) # Save frequency analysis results to hdf5
        if alltime: # Since we're plotting over all time, no need to normalize to anything else. Just make the plots.
            plotme(data2, outdir=outdir, pltdict=pltdict, fld_id = fld_id, alltime=alltime) # Make plots and save to png
        data2_dict[fld_id] = data2
        pltdict_dict[fld_id] = pltdict

    return data2_dict, pltdict_dict

def freqanalyze(fns, data=None, datpath=None, fld_id = 'Ez', divsp = 1, pool = None):
    # Inputs:
    #   divsp: integer, divisor by which to reduce the spatial resolution (e.g. divsp = 2 reduces field dimensions from 300x200 to 150x100)
    #   data: (Optional) A python dictionary returned by h5stitch2D.py "fields2D()". If none is supplied, it will be read in by reading from datah5 (fields2D.hdf5) or by looking at p4 filenames.
    #   datpath: (Optional) A string path of the hdf5 file from which data can be loaded.
    #   fns: list of filenames to analyze, if neither data nor datah5 are given.
    #   pool: The pool multiprocessing threads to use when doing the FFT step. If pool = None (default), just use serial processing.
    #   fld_id: A string specifying the field component to analyze. This field must be contained within the "data" dictionary

    # Decide how to read in the data
    if data:
        #Data supplied at input to freqanalyze(), so no need to read in the files.
        print "Extracting fields data directly from python Data dict."
        times = data['times']*1e6 # times, converted to fs
        fns = data['filenames']
        xgv = data['xgv'][::divsp]*1e4 # spatial X, converted to microns
        zgv = data['zgv'][::divsp]*1e4 # spatial Z, converted to microns
        Ez = data[fld_id][:,::divsp,::divsp] ## I call it "Ez" as a variable name, but this could be any field.
    elif datpath:
        print "Reading fields2D data from the HDF5 file."
        with h5py.File(datpath, 'r') as data:
            times = data['times'][...]*1e6 # times, converted to fs
            fns = data['filenames'][...]
            xgv = data['xgv'][::divsp][...]*1e4 # spatial X, converted to microns
            zgv = data['zgv'][::divsp][...]*1e4 # spatial Z, converted to microns
            Ez = data[fld_id][:,::divsp,::divsp][...] ## I call it "Ez" as a variable name, but this could be any field.
    else:
        print "Reading ", len(fns), "files."
        # Load in the data
        data = fields2D(fns, fld_ids = [fld_id], pool = pool)
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

    # Assume equal spacing in time and space, and get the deltas
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

    #print "Doing calculations"
    if pool:    
        print "Parallel FFT started."
        Eft = fftpar(Ez, pool) # PARALLEL FFT OPTION
    else:
        print "Serial FFT started."
        Eft = np.fft.rfft(Ez, axis = 0) # SERIAL FFT OPTION

    # Output a variety of analyzed data for plotting
    data2['Imap'] = np.sum(Ez**2,0)
    pwr = np.absolute(Eft)**2
    data2['FTmap_0_5'] = np.sum(pwr[np.logical_and(freq > 0.3, freq < 0.7)],0)
    data2['FTmap_1'] = np.sum(pwr[np.logical_and(freq > 0.9, freq < 1.1)],0)
    data2['FTmap_1_5'] = np.sum(pwr[np.logical_and(freq > 1.3, freq < 1.7)],0)
    data2['FTmap_2'] = np.sum(pwr[np.logical_and(freq > 1.8, freq < 2.3)],0)
    data2['FTmap_0_85'] = np.sum(pwr[np.logical_and(freq > 0.7, freq < 0.9)],0)
    data2['pwr_sum'] = np.sum(pwr, (1,2))

    return data2

def myfig(data2,mapID,pltdict,xgv,zgv,tstring,fld_id,sticker,title,fignum):
    """A very custom plot for color plots of field values cut at frequency bands."""
    fig = plt.figure(fignum)
    plt.clf() # Clear the figure
    C = data2[mapID]
    if pltdict:
        cmin = pltdict[mapID]['min']
        cmax = pltdict[mapID]['max']
    else:
        cmin = np.min(C)
        cmax = np.max(C)
    ax = plt.subplot(111)
    im = ax.pcolorfast(xgv,zgv,C, cmap='viridis')
    ax.set_xlabel(r'X ($\mu m$)')
    ax.set_ylabel(r'Z ($\mu m$)')
    ax.set_title(title, fontsize=20)
    ax.text(np.min(xgv) + 1, np.max(zgv) - 3, tstring, fontsize=24, color='white')
    ax.text(np.max(xgv) - 6, np.min(zgv) + 2, fld_id, fontsize=44, color='white')
    ax.text(np.min(xgv) + 1, np.min(zgv) + 2, sticker, fontsize=44, color='white')
    cbar = fig.colorbar(im, label=r'Power density (a.u.)')
    im.set_clim(vmin=cmin, vmax=cmax)
    return fig

def plotme(data2, outdir='', pltdict = None, fld_id = 'Fld', alltime=False):
    """ Plots a variety of parameters found in the outputs of frequency analysis, then saves a png of each."""
    xgv = data2['xgv_um']
    zgv = data2['zgv_um']
    freq = data2['freq']
    times = data2['times_fs']
    pwr_sum = data2['pwr_sum']

    if alltime: # If this is the all-time rather than time-resolved analysis, put a different label on plot
        maxtime = np.max(times)
        mintime = np.min(times)
        tstring = "All time"
    else: # Otherwise, put the standard "t = XX fs +/- 20 fs" label.
        meantime = np.mean(times)
        tplus = np.max(times) - np.mean(times)
        tstring = 't=' + "{:.1f}".format(meantime) + " fs $\pm$ " + "{:.1f}".format(tplus) + " fs"

    # Make some plots
    #print "Making plots."
    figs = [] # Figure handles
    labs = [] # Plot name labels

    # Power spectrum
    label = 'Power spectrum'
    fig = plt.figure(0)
    figs.append(fig) # Add the figure to the figures list
    labs.append(label)
    plt.clf() # Clear the figure
    if pltdict:
        ymax = np.log10(pltdict['pwr_sum']['max'])
        ymin = ymax - 3 # Give 3 orders of magnitude
    else:
        ymax = np.log10(np.max(pwr_sum))
        ymin = ymax - 3 # Give 3 orders of magnitude
    xmin = 0
    xmax = 2.5
    ax = plt.subplot(111)
    ax.plot(freq, np.log10(pwr_sum))
    ax.set_xlabel('Frequency (normalized to laser fundamental)')
    ax.set_ylabel('Log$_{10}$(Spectral power density) (a.u.)')
    ax.set_title('Power spectrum (Log scale), ' + tstring, fontsize=16)
    ax.set_xlim(xmin, xmax)
    ax.xaxis.grid() # vertical lines
    ax.text(xmax - (xmax - xmin)/7, ymin + (ymax - ymin)/16, fld_id, fontsize=44, color='black')
    ax.set_ylim(ymin, ymax)

    # All frequencies map
    label = 'Intensity'
    sticker = r'$I$'
    title = fld_id + r', Intensity map'
    mapID = 'Imap'

    fig = myfig(data2,mapID,pltdict,xgv,zgv,tstring,fld_id,sticker,title,1)
    figs.append(fig)
    labs.append(label)

    
    # Half omega map
    label = 'Half omega'
    sticker = r'$\omega/2$'
    title = fld_id + r', 0.3 to 0.7 $\omega$'
    mapID = 'FTmap_0_5'

    fig = myfig(data2,mapID,pltdict,xgv,zgv,tstring,fld_id,sticker,title,2)
    figs.append(fig)
    labs.append(label)


    # Omega map
    label = 'Omega'
    sticker = r'$1\omega$'
    title = fld_id + r', 0.9 to 1.1 $\omega$'
    mapID = 'FTmap_1'

    fig = myfig(data2,mapID,pltdict,xgv,zgv,tstring,fld_id,sticker,title,3)
    figs.append(fig)
    labs.append(label)


    # Three-halves omega map
    label = 'Three-halves omega'
    sticker = r'$3\omega/2$'
    title = fld_id + r', 1.3 to 1.7 $\omega$'
    mapID = 'FTmap_1_5'

    fig = myfig(data2,mapID,pltdict,xgv,zgv,tstring,fld_id,sticker,title,4)
    figs.append(fig)
    labs.append(label)

    # Two omega map
    label = 'Two omega'
    sticker = r'$2\omega$'
    title = fld_id + r', 1.8 to 2.3 $\omega$'
    mapID = 'FTmap_2'

    fig = myfig(data2,mapID,pltdict,xgv,zgv,tstring,fld_id,sticker,title,5)
    figs.append(fig)
    labs.append(label)


    # 0.85x omega map
    label = '0_85 omega'
    sticker = r'$0.85\omega$'
    title = fld_id + r', 0.7 to 0.9 $\omega$'
    mapID = 'FTmap_0_85'

    fig = myfig(data2,mapID,pltdict,xgv,zgv,tstring,fld_id,sticker,title,6)
    figs.append(fig)
    labs.append(label)


    print "Saving figures"
    # Save the figures into custom folders and with time-based filenames.
    if alltime: # If this is the "All times" analysis, label the file with "00000.*"
        tlabel = ''.zfill(5)
    else: # Otherwise, do the normal labeling
        tlabel = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the file label be '00512.*' for t = 51.2342 fs

    plotdir = subdir(outdir, 'FreqPlots')
    flddir = subdir(plotdir, fld_id)
    for i in range(len(figs)):
        figdir = subdir(flddir, fld_id + ' -' + str(i) + '- ' + labs[i])
        fn = os.path.join(figdir, tlabel + '.png') # For example, store the figure as '[outdir]/FreqPlots/Ex/Ex -2- Half omega/00512.png'
        figs[i].savefig(fn)

    #plt.close('all') # This is nice but not necessary, since I numbered the figures (they will be re-used across calls)

def fftpar(Ei, pool):
    """ Parallel temporal fast fourier transform on a (time x space x space) field array. Returns a (freq x space x space) array."""
    # It may not actually be all that efficient to spin off a bunch of threads for this calculation, but it is still a lot faster than one thread.
    # Prior to calling, get your pool via: pool = Pool(10) for a 10-thread pool
    return np.swapaxes(np.array((pool.map(np.fft.rfft, Ei.swapaxes(0,2)))),0,2)

def subdir(folder, name):
    """ Make a subdirectory in the specified folder, if it doesn't already exist"""
    subpath = os.path.join(folder,name)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    return subpath

def freqSave(data2, outdir = '', fld_id = '', alltime=False):
    """Save the data used in frequency plots to an hdf5 file."""
    if alltime: # If this is the "All times" analysis, label the file with "00000.*"
        tlabel = ''.zfill(5)
    else: # Otherwise, do the normal labeling
        meantime = np.mean(data2['times_fs'])
        tlabel = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the file label be '00512.*' for t = 51.2342 fs

    plotdir = subdir(outdir, 'FreqData')
    flddir = subdir(plotdir, fld_id)
    h5path = saveData(data2, flddir, tlabel + '.hdf5') # E.g. save the file as '[outdir]/FreqData/Ex/00512.hdf5'
    return h5path

    
def saveData(data2, folder = '', name = 'data2.hdf5'):
    """Save a python dictionary containing only top-level NumPy arrays (example: data2 in this document) to an HDF5 file, with moderate gzip compression"""
    print "Saving HDF5 file."
    h5path = os.path.join(folder, name)
    with h5py.File(h5path, 'w') as f:
        for k in data2:
            f.create_dataset(k, data = data2[k], compression='gzip', compression_opts=4)

    return h5path

def loadData(h5path):
    """Read a simple HDF5 file containing only top-level arrays back into a python dict of NumPy arrays"""
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

