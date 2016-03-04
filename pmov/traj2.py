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
import time
from multiprocessing import Pool

try:
    from mpi4py import MPI
except:
    print "WARNING: MPI4PY FAILED TO LOAD. DO NOT CALL PARALLEL MPI FUNCTIONS."


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
    
def mpiTraj(p4dir, h5fn = None):
    """ Assume we have greater than one processor. Rank 0 will do the hdf5 stuff"""
    # Set some basic MPI variables
    nprocs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()
    comm = MPI.COMM_WORLD

    # Check that we have at least two processors in our MPI setup    
    if nprocs < 2:
        raise Exception("Not enough processors. Need at least two, or call some other function.")

    # Everyone, get your bearings on the task to be performed
    # Get a sorted list of filenames for the pmovie files
    fns = ls.getp4(p4dir, prefix = "pmovie") # Get list of pmovieXX.p4 files, and get the list sorted in ascending time step order
    nframes = len(fns) # Number of pmovie frames to prepare for

    # Rank0: Read in the first file and get reference data, spread that around  
    if rank == 0: # Rank 0, start working on that HDF5
        # Read the first frame, to get info like the length
        print "Rank 0 speaking, I'm reading in the first frame. Everyone else sit tight."
        data_ref, _ = sortOne(fns[0])
        print "Ok, I'm going to spread that around, now."
    else:
        data_ref = None
    data_ref = comm.bcast(data_ref, root=0)
    
    nparts = len(data_ref) # Number of particles extracted from this first frame (which will determine the rest, as well)

    # Assign files to each processor
    framechunks = sf.chunkFair(range(nframes), (nprocs - 1))
    myframes = np.array(framechunks[rank - 1], dtype='i') #TODO: Better ordering of chunks (non-sequential) # If there are 4 processors, break into 3 chunks (as 0th processor just writes files)
    
    if rank == 0: # Rank 0, start working on that HDF5
        print "Good stuff. I'm going to get started on this HDF5; I'll let you know how it goes."
        # Determine the filename for output of hdf5 trajectories file
        if not h5fn:    
            h5fn = os.path.join(p4dir, 'traj.h5')

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
                print "Rank 0: Collecting file ", i, " of ", nframes   
                datdict = comm.recv(source=MPI.ANY_SOURCE, tag=i)  # Retrieve the result
                datnew = datdict['datnew']
                stats = datdict['stats']
                goodcdt = datdict['goodcdt']
                badcdt = np.logical_not(goodcdt) # Flip the sign of good condit
                datnew[badcdt] = data_ref[badcdt] # Fill in the missing particles
                f['t'][i] = stats['t']
                f['step'][i] = stats['step']
                f['gone'][i] = badcdt # Flag the particles that were missing
                for k in goodkeys:
                    f[k][i] = datnew[k]
                
                data_ref = datnew # The new array becomes the reference for next iteration
            print "ELAPSED TIME (secs)", (dt.now() - t1).total_seconds()

    else:
        print "I am a servant of the ranks."
        print "Rank", rank, ": I have frames:", myframes
        for i in range(len(myframes)):
            ix = myframes[i] # The index that refers to the entire list (of all files)
            print "Rank", rank, ": working on file", ix
            fn = fns[ix]
            datnew, stats, goodcdt = crunchOne(fn, data_ref)
            datdict = {}
            datdict['datnew'] = datnew
            datdict['stats'] = stats
            datdict['goodcdt'] = goodcdt
            comm.send(datdict, dest=0, tag=ix) # Note: comm.isend would give an EOFError, for some reason, so don't use it.
