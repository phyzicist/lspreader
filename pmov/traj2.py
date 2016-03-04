# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 17:31:13 2016

Scott's implementation of pmovie sorting

@author: Scott
"""


import lspreader2 as rd
import sftools as sf
import lstools as ls
import numpy as np
import h5py
import os
from datetime import datetime as dt

from multiprocessing import Pool


def sortOne(fn, npz=False):
    """ Sorts one pmovie file according to initial position (first x, then z, then y).
    Assumes there is only one time step per pmovie file. Otherwise, takes only the first of these time steps.
    Inputs:
        fn: string, filename e.g.  fn="../../pmovie0004.p4(.gz)"
        npz: bool, True if we would like to save the sorted array as pmovieXX.p4(.gz).npy, in the same folder as the original file
    Outputs:
        data: 1D numpy array with several dtype fields, equivalent to the output of one frame of lspreader's read_pmovie(), but sorted by initial particle position
        stats: A dict containing information such as step number, step time, and number of particles
        (Optional) Saves the above array to file (off by default)        
    """
    
    frames = rd.read_movie2(fn) # Read in the pmovie file using lspreader2
    data = frames[0]['data'] # Extract the first (and assumed only) time step of particles from the pmovie file
    
    del frames[0]['data'] # Remove the reference to the data array from the frame dictionary
    stats = frames[0] # Refer to this array as the "stats" dictionary, which no longer contains the data
    
    data = np.sort(data, order=['xi','zi','yi']) # Sort the array in order of particle initial position: first x, then z, then y
    
    if npz:
        np.savez(fn + ".npz", data, np.array(stats)) # If the npy flag was set to True at input, save this to as a numpy array to file     
    
    return data, stats # Return the sorted array

def fillGaps(data, data_ref):
    """ Fill gaps (missing particles) in "data" with particles from "data_ref". Fill in missing particles in data. That is, those present in data_ref but not in data will be taken from data_ref and inserted into data.
    Inputs:
        Note - inputs and outputs are all 1D numpy array with several dtype fields, equivalent to the output of one frame of lspreader's read_pmovie(), but sorted by initial particle position
        data: sorted data array with missing particles
        data_ref: sorted data array with all particles (the template, e.g. from previous time step) (equal to in length or longer than data)
    Outputs:
        data_new: sorted data array which is a blend of data_ref and data, where missing particles have been replaced.
    """
    goodcdt = np.zeros(data_ref.shape, dtype=bool) # Allocate an array where present (not missing) particles will be marked as True. Initialize to False (all particles missing).
    ix1 = 0 # The index for data_ref
    ix2 = 0 # The index for data (will always be equal to or less than ix1)
    while ix2 < len(data):
        # Iterate over each of the reference particles, checking if they are missing from the data.
        # For each reference particle, check if any xi, zi, or yi are unequal
        if (data[ix2]['xi'] != data_ref[ix1]['xi']) & (data[ix2]['zi'] != data_ref[ix1]['zi']) & (data[ix2]['yi'] != data_ref[ix1]['yi']): # If xi, yi, zi aren't equal, this is a missing particle
            pass
        else: # This is not a missing particle.
            goodcdt[ix1] = True # Mark this particle as present
            ix2 += 1 # Move onto the next data particle only if the particle wasn't missing. Otherwise, stay on this data particle.
        ix1 += 1 # Every iteration, move on to the next reference particle
    
    # Copy into the new array, dat3s. Could also just do a rolling update.
    data_new = np.copy(data_ref) # Make a deep copy of the reference array, as your new output array
    data_new[goodcdt] = data # Fill in all particles that were present in your data (the missing particles will stay what they were in the reference data)
    
    print "Fraction missing:", 1 - float(len(goodcdt[goodcdt]))/float(len(goodcdt))

    return data_new, goodcdt  # Return the data array, with all particles now present (gaps are filled)

def crunchOne(fn, data_ref):
    """ Open file, fill in missing particles (according to data_ref). Consolidated into a single function for use in parallel processing """
    dattmp, stats = sortOne(fn)
    data_new, goodcdt = fillGaps(dattmp, data_ref)
    return data_new, stats, goodcdt

def serTraj(p4dir, h5fn = None):
    """ Get the trajectories of a folder into an HDF5 file (in serial)
    Inputs:
        p4dir: string, path to the containing p4movie files
        h5fn: (Optional) string, the full-path filename (e.g. h5fn = ".../.../mytrajs.h5"
    """
    
    # Determine the filename for output of hdf5 trajectories file
    if not h5fn:    
        h5fn = os.path.join(p4dir, 'traj.h5')

    # Get a sorted list of filenames for the pmovie files
    fns = ls.getp4(p4dir, prefix = "pmovie") # Get list of pmovieXX.p4 files, and get the list sorted in ascending time step order
        
    nframes = len(fns) # Number of pmovie frames to prepare for
    # Read the first frame, to get info like the length
    datref, _ = sortOne(fns[0])
    
    nparts = len(datref) # Number of particles extracted from this first frame (which will determine the rest, as well)
    
    t1 = dt.now()
    goodkeys = ['xi', 'zi', 'x', 'z', 'ux', 'uy', 'uz']
    # Open the HDF5 file, and step over the p4 files
    with h5py.File(h5fn, "w") as f:
        # Allocate the HDF5 datasets
        f.create_dataset("t", (nframes,), dtype='f')
        f.create_dataset("step", (nframes,), dtype='i')
        for k in goodkeys:
            f.create_dataset(k, (nframes, nparts,), dtype='f')
        
        # Now, get into the beef of stepping over the files
        for i in range(nframes):
            print "Reading file ", i, " of ", nframes
            dattmp, stats = sortOne(fns[i])
            
            datnew, goodcdt = fillGaps(dattmp, datref)
            
            f['t'][i] = stats['t']
            f['step'][i] = stats['step']
            for k in goodkeys:
                f[k][i] = datnew[k]
                
            datref = datnew # The new array becomes the reference for next iteration
    t2 = dt.now()
    delta = t1 - t2
    print "Elapsed time", delta.total_seconds()
    return datref

def parTraj(p4dir, h5fn = None, nprocs = 4):
    """ Get the trajectories of a folder into an HDF5 file (in parallel)
    Inputs:
        p4dir: string, path to the containing p4movie files
        h5fn: (Optional) string, the full-path filename (e.g. h5fn = ".../.../mytrajs.h5"
        nprocs: The number of processors to open up with the multiprocessing unit Pool
    """
    
    # Determine the filename for output of hdf5 trajectories file
    if not h5fn:    
        h5fn = os.path.join(p4dir, 'traj.h5')

    pool = Pool(processes=nprocs)    # start nprocs worker processes

    # Get a sorted list of filenames for the pmovie files
    fns = ls.getp4(p4dir, prefix = "pmovie") # Get list of pmovieXX.p4 files, and get the list sorted in ascending time step order
        
    nframes = len(fns) # Number of pmovie frames to prepare for
    # Read the first frame, to get info like the length
    data_ref, _ = sortOne(fns[0])
    
    nparts = len(data_ref) # Number of particles extracted from this first frame (which will determine the rest, as well)


    # Now, get into the beef of stepping over the files
    print "First pass: Assigning multiple tasks."
    tasks = [None]*nframes # List of assigned processor tasks
    for i in range(nframes):
        tasks[i] = pool.apply_async(crunchOne, (fns[i], data_ref))
    
    t1 = dt.now()
    goodkeys = ['xi', 'zi', 'x', 'z', 'ux', 'uy', 'uz']
    # Open the HDF5 file, and step over the p4 files
    with h5py.File(h5fn, "w") as f:
        # Allocate the HDF5 datasets
        f.create_dataset("t", (nframes,), dtype='f', compression="gzip")
        f.create_dataset("step", (nframes,), dtype='int32', compression="gzip")
        f.create_dataset("gone", (nframes, nparts,), dtype='bool', compression="gzip")
        for k in goodkeys:
            f.create_dataset(k, (nframes, nparts,), dtype='f', compression="gzip")
        
        # Now, iterate over the files (collect their data and save to the HDF5)
        for i in range(nframes):
            print "Second pass: Collecting file ", i, " of ", nframes
            
            datnew, stats, goodcdt = tasks[i].get(timeout=60*60) # Retrieve the result
            badcdt = np.logical_not(goodcdt) # Flip the sign of good condit
            datnew[badcdt] = data_ref[badcdt] # Fill in the missing particles
            f['t'][i] = stats['t']
            f['step'][i] = stats['step']
            f['gone'][i] = badcdt # Flag the particles that were missing
            for k in goodkeys:
                f[k][i] = datnew[k]
            
            data_ref = datnew # The new array becomes the reference for next iteration
    t2 = dt.now()
    delta = t1 - t2
    print "Elapsed time", delta.total_seconds()
    pool.close() # Close the parallel pool
    return data_ref