# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:19:17 2016

Functions no longer in use. Retained just in case.

@author: Scott
"""

############ H5STITCH2D FUNCTIONS

import h5py
import os
import lspreader2 as rd
import numpy as np
import re
import gzip
from multiprocessing import Pool
import sys
from sftools import chunkIt

try:
    from mpi4py import MPI
except:
    print "WARNING: MPI4PY FAILED TO LOAD. DO NOT CALL PARALLEL FUNCTIONS."

## DINOSAUR FUNCTIONS!!! (MPI-BASED) DO NOT USE:
def h5fields2Dser(folder, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz']):
    """ 100% serial processing. Read in a folder full of flds*.p4(.gz) files, stitching them together assuming 2D assumptions, and create an HDF5 file """
    
    # Extract "E" from "Ex" in fld_ids: fields (for read_flds2() call later)    
    flds = list(set(s[:-1] for s in fld_ids)) # we strip the "x","y","z" last character off our fields, then set() gives only unique elements of a list, and list() converts this set back to list
    # E.g. if fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], flds = ['E','B']
    
    if not h5path:
        h5path = os.path.join(folder, 'fields2D.hdf5')

    fns = getfnsp4(folder, 'flds')
    nfiles = len(fns)

    
    # Open the HDF5 file
    with h5py.File(h5path,'w') as f:
        # Read in the first file as a template
        doms, header = rd.read_flds2(fns[0], flds=flds)
        
        # Pre-allocate the HDF5 arrays
        _, xgv, zgv = rd.stitch2D(doms, fld_ids[0]) # Extract the interesting fields and stitch together for a template
       
        f.create_dataset('times', (nfiles,), compression='gzip', compression_opts=4, dtype='float32')
        f.create_dataset('xgv', compression='gzip', compression_opts=4, data = xgv)
        f.create_dataset('zgv', compression='gzip', compression_opts=4, data = zgv)
        f.create_dataset('filenames', data = fns) # A bit overkill to store the full filenames, but who cares? Space is not too bad.
        for k in fld_ids:
            f.create_dataset(k,(nfiles,len(zgv),len(xgv)), compression='gzip', compression_opts=4, dtype='float32')
        
        # Read in all the data and fill up the HDF5 file
        for i in range(nfiles): # Iterate over the files
            fn = fns[i]
            doms, header = rd.read_flds2(fn, flds=flds)
            f['times'][i] = header['timestamp']
            for fld_id in fld_ids: # Iterate over the requested fields, stitching then adding them to HDF5 file
                fld, xgv, zgv = rd.stitch2D(doms, fld_id)
                f[fld_id][i,:,:] = fld
    return h5path


def h5fields2Dpar(folder, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], flds = ['E', 'B'], compress=True):
    """ Parallel processing, with serial HDF5. Read in a folder full of flds*.p4 file, stitching them together assuming 2D assumptions, and create an HDF5 file """
    # Note that compression filters do not work with parallel I/O as of h5py version 2.5.0. Hence, an extra obnoxious step at the end with a single processor re-saving the files (unless compress=False flag set, which would make the computer time faster)
    if not h5path:
        h5path = os.path.join(folder, 'fields2D.hdf5')

    fns = getfnsp4(folder, 'flds')
    nfiles = len(fns)

    nprocs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()

    comm = MPI.COMM_WORLD

    if rank == 0:
        idx_all = range(nfiles) # This list of indices can be split up later across processors, outside this function
        idx_chunked = chunkIt(idx_all, nprocs) # Break the list into N chunks, where N = number of processors
        print "Scattering the fld*.p4 files across " + str(nprocs) + " nodes."
    else:
        idx_chunked = None
    
    idx_part = comm.scatter(idx_chunked,root=0) # Assign each processor one chunk of the list. Therefore, idx_part (list of indices) specifies which elements to analyze here.


    # Read in all the data and create temporary lists, 'times' and 'flds'
    times_part = []
    flds_part = []       
    for i in idx_part: # Iterate over the files for this processor
        fn = fns[i]
        doms, header = rd.read_flds2(fn, flds=flds)
        times_part.append(header['timestamp'])
        fld = {}
        for fld_id in fld_ids: # Iterate over the requested fields, stitching then adding them to HDF5 file
            fld[fld_id], xgv, zgv = rd.stitch2D(doms, fld_id)
        flds_part.append(fld)
    #print "Succeeded in adding " + str(len(idx_part)) + " files on processor " + str(rank)

    flds_complete = comm.gather(flds_part, root=0)
    idx_complete = comm.gather(idx_part, root=0)
    times_complete = comm.gather(times_part, root=0)
    if rank == 0:
        print "Data read-in complete, and data transferred to proc. 0. Stitching into HDF5 on proc. 0."
        # Use the 0th core to copy everything into a compressed dataset.  

        # Read in the first file as a template
        doms, header = rd.read_flds2(fns[0], flds=flds)
        
        # Get xgv and zgv
        _, xgv, zgv = rd.stitch2D(doms, fld_ids[0]) # Extract the interesting fields and stitch together for a template
       

        # Consolidate data into appropriate times and fields arrays
        times = np.zeros((nfiles,), dtype='float32')
        fld_array = {}     
        for fld_id in fld_ids:        
            fld_array[fld_id] = np.zeros((nfiles,len(zgv),len(xgv)), dtype='float32')

        for i in range(nprocs):
            for j in range(len(idx_complete[i])):
                idx = idx_complete[i][j]
                times[idx] = times_complete[i][j]
                for fld_id in fld_ids:
                    fld_array[fld_id][idx,:,:] = flds_complete[i][j][fld_id]

        # Open the HDF5 file and write the data
        with h5py.File(h5path,'w') as f:
            f.create_dataset('xgv', data = xgv, compression='gzip')
            f.create_dataset('zgv', data = zgv, compression='gzip')
            f.create_dataset('filenames', data = fns, compression='gzip') # A bit overkill to store the full filenames, but who cares? Space is not too bad.
            f_times = f.create_dataset('times', data = times, compression='gzip')

            for k in fld_ids:
                f.create_dataset(k, data = fld_array[k], compression='gzip')

    return h5path


def h5fields2Dpar2(folder, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz']):
    """ Parallel processing. DOES NOT WORK IN CURRENT FORM. Creates the HDF5 file with less pain. However, you need lots of RAM on node0 to hold all the data."""
    if not h5path:
        h5path = os.path.join(folder, 'fields2D.hdf5')
    fns = getfnsp4(folder)
    nfiles = len(fns)

    # Do some fancy parallel stuff
    nprocs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()

    comm = MPI.COMM_WORLD

    # Split the files across nodes
    if rank == 0:
        fns_chunked = chunkIt(fns, nprocs) # Break the list into N chunks, where N = number of processors
        print "Scattering " + str(nfiles) + " fld*.p4 files across " + str(nprocs) + " nodes, to read into NumPy arrays and pass back."
    else:
        fns_chunked = None
    
    fns_part = comm.scatter(fns_chunked,root=0) # Assign each processor one chunk of the files. Therefore, each processor can look at fns_part (partial list of filenames) as a list of files for it to analyze.

    # Read files assigned to this processor into RAM
    data_part = fields2D(fns_part, fld_ids=fld_ids)
    # Perform any additional analysis

    ## Consolidate data on the 0th core to copy everything into a compressed dataset.      
    data_chunked = comm.gather(data_part, root=0)

    
    print data_chunked[0]
    if rank == 0:
        print "Data read-in complete, and data transferred to proc. 0. Stitching into HDF5 on proc. 0."
        
        with h5py.File(h5path,'w') as f:
            print "HDF5 file opened."
            # Build the HDF5 file, assuming every element in "data" is a NumPy array
            print "Saving arrays in RAM to HDF5"

            # Choose a non-empty chunk to extract keys
            count = 0
            while (len(fns_chunked[count]) < 1):
               print 'The count is:', count
               count = count + 1
            print "Chose test idx of", count

            data_test = data_chunked[count]
            # Pre-allocate HDF5 arrays
            nostretch = ['xgv', 'ygv', 'zgv'] # keys of arrays that should not scale with the filenumber (array built from one file)
            stretch = [] # keys of arrays that will scale with the filenumber (array built from many files); built in next loop
            for k in data_test:
                shape = data_test[k].shape
                if k in nostretch:
                    #print 'Will not stretch to number of files: ', k
                    shape_new = shape
                else:
                    stretch.append(k)
                    #print 'Will stretch to number of files:', k
                    shape_list = list(shape)
                    shape_list[0] = nfiles # resize the test shape to account for all files
                    shape_new = tuple(shape_list)
                print k, ":", shape, " resized=> ", shape_new
                f.create_dataset(k, shape_new, dtype=data_test[k].dtype, compression='gzip', compression_opts=4)

            # Fill in the HDF5 arrays
            nsaved = 0
            for i in range(nprocs):
                nfiles = len(fns_chunked[i])
                print nfiles, "files and "
                for j in range(nfiles):
                    print 'J', j, 'J + nsaved', j + nsaved
                    for k in stretch:
                        print k
                        print data_chunked[i][k].shape, f[k].shape
                        if k == 'times':
                            print 'Time', data_chunked[i][k][j], '===>', f[k][j + nsaved]
                            print k, data_chunked[i][k].shape, '===>', f[k].shape
                        else:
                            print k, data_chunked[i][k].shape, '===>', f[k].shape
                        f[k][j + nsaved] = data_chunked[i][k][j]
                        print k, data_chunked[i][k].shape, '===>', f[k].shape
                    print "Saved one whole file!"
                nsaved = nsaved + nfiles
                print "Finished data for processor", i, "Have now saved", nsaved, "files."

            print "All done!"



########### FREQUENCY ANALYSIS FUNCTIONS
import h5py
import os
import numpy as np
from h5stitch2D import getfnsp4, fields2D
from sftools import subdir, chunkFx, chunkIt

# Matplotlib stuff
import matplotlib as mpl
#mpl.use('Agg') # Let's call this at an earlier point, shall we?
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

## DINOSAUR FUNCTIONS NO LONGER IN USE


def freqanalyze(fns, data=None, datpath=None, fld_id = 'Ez', divsp = 1, pool = None):
    """ Analyze a list of filenames (alternatively, already loaded data) and output frequency plots, etc. to file.
    Inputs:
        divsp: integer, divisor by which to reduce the spatial resolution (e.g. divsp = 2 reduces field dimensions from 300x200 to 150x100)
        data: (Optional) A python dictionary returned by lstools.py "ls.fields2D()". If none is supplied, it will be read in by reading from datah5 (fields2D.hdf5) or by looking at p4 filenames.
        datpath: (Optional) A string path of the hdf5 file from which data can be loaded.
        fns: list of filenames to analyze, if neither data nor datah5 are given.
        pool: The pool multiprocessing threads to use when doing the FFT step. If pool = None (default), just use serial processing.
        fld_id: A string specifying the field component to analyze. This field must be contained within the "data" dictionary
    """
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
        data = ls.fields2D(fns, fld_ids = [fld_id], pool = pool)
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
def mainPar():
    """DINOSAUR FUNCTION, no longer in use. Retained as an example of using MPI."""
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
    """DINOSAUR function, no longer in use. Retained as an example of not using MPI."""
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
