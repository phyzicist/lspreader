# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 17:31:13 2016

Scott's implementation of pmovie sorting

@author: Scott

Changelog:
2016-07-06 Updated hashing to use Gregory's. Also, now fills in particles with their initial conditions (rather than most recent).

TODO:
* Compute kinetic energy and angle for each particle at each step and add to the HDF5

"""


import lspreader2 as rd
import sftools as sf
import lstools as ls
import numpy as np
import h5py
import os
from datetime import datetime as dt
import time
import scipy.constants as sc
import traceback
import numpy.lib.recfunctions as rfn
try:
    from mpi4py import MPI
except:
    print "WARNING: MPI4PY FAILED TO LOAD. DO NOT CALL PARALLEL MPI FUNCTIONS."

def firsthash(frame,dims,removedupes=False):
    '''
    [Copied from Gregory Ngirmang's lspreader.
    https://github.com/noobermin/lspreader]
    Hashes the first time step. Only will work as long as
    the hash can fit in a i8.
    Parameters:
    -----------
      frame : first frame.
      dims  :  iterable of strings for dimensions.
    Keywords:
    ---------
      removedupes: specify duplicates for the given frame.
    
    Returns a dictionary of everything needed
    to generate hashes from the genhash function.
    
    '''
    #hashes must have i8 available
    #overwise, we'll have overflow
    def avgdiff(d):
        d=np.sort(d);
        d = d[1:] - d[:-1]
        return np.average(d[np.nonzero(d)]);
    ip    = np.array([frame['data'][l] for l in dims]).T;
    avgdiffs = np.array([avgdiff(a) for a in ip.T]);
    mins  = ip.min(axis=0);
    ips = (((ip - mins)/avgdiffs).round().astype('i8'))
    pws  = np.floor(np.log10(ips.max(axis=0))).astype('i8')+1
    pws = list(pws);
    pw = [0]+[ ipw+jpw for ipw,jpw in
               zip([0]+pws[:-1],pws[:-1]) ];
    pw = 10**np.array(pw);
    #the dictionary used for hashing
    d=dict(labels=dims, mins=mins, avgdiffs=avgdiffs, pw=pw);
    if removedupes:
        hashes = genhash(frame,d,removedupes=False);
        #consider if the negation of this is faster for genhash
        uni,counts = np.unique(hashes,return_counts=True);
        d.update({'dupes': uni[counts>1]})
    return d;

def genhash(frame,d,removedupes=False):
    '''
    [Copied from Gregory Ngirmang's lspreader.
    https://github.com/noobermin/lspreader]
    Generate the hashes for the given frame for a specification
    given in the dictionary d returned from firsthash.
    Parameters:
    -----------
      frame :  frame to hash.
      d     :  hash specification generated from firsthash.
    Keywords:
    ---------
      removedupes: put -1 in duplicates
    
    Returns an array of the shape of the frames with hashes.
    '''
    ip = np.array([frame['data'][l] for l in d['labels']]).T;
    scaled = ((ip - d['mins'])/d['avgdiffs']).round().astype('i8');
    hashes = (scaled*d['pw']).sum(axis=1);
    #marking duplicated particles
    if removedupes:
        dups = np.in1d(hashes,d['dupes'])
        hashes[dups] = -1
    return hashes;
    
def addhash(frame,d,removedupes=False):
    '''
    [Copied from Gregory Ngirmang's lspreader.
    https://github.com/noobermin/lspreader]
    helper function to add hashes to the given frame
    given in the dictionary d returned from firsthash.
    Parameters:
    -----------
      frame :  frame to hash.
      d     :  hash specification generated from firsthash.
    Keywords:
    ---------
      removedupes: put -1 in duplicates
    
    Returns frame with added hashes, although it will be added in
    place.
    '''
    hashes = genhash(frame,d,removedupes);
    frame['data'] = rfn.rec_append_fields(
        frame['data'],'hash',hashes);
    return frame;

def sortOne(fn, hashd=None):
    """ Loads one pmovie file and hashes it, then sorts it by hash. Duplicately-hashed particles are deleted from the output list, with prejudice.
    Assumes there is only one time step per pmovie file. Otherwise, takes only the first of these time steps.
    Inputs:
        fn: string, filename e.g.  fn="../../pmovie0004.p4(.gz)"
        hashd: The hashing dictionary greated by genhash(). If None, assume this is the first frame and defines future hashing.
    Outputs:
        data: 1D numpy array with several dtype fields, equivalent to the output of one frame of lspreader's read_pmovie(), but sorted by initial particle position
        stats: A dict containing information such as step number, step time, and number of particles
    """
    frame = rd.read_movie2(fn)[0] # Read in the pmovie file using lspreader2. Assume first frame is the only frame.
    print len(frame['data'])
    if hashd is None: # First pmovie file; define the hashing functions/parameters into a dict called "hashd"
        hashd = firsthash(frame, ['xi','zi','yi'], removedupes=True)

    # Splice "hash" field as a frame['data'] field
    frame = addhash(frame, hashd, removedupes=True)

    # Sort by hash
    frame['data'].sort(order='hash')

    data = frame['data'][frame['data']['hash'] != -1]
    print "Fraction of duplicates:", (0.0 + len(frame['data']) - len(data))/len(frame['data'])

    del frame['data'] # Remove the reference to the data array from the frame dictionary
    stats = frame # Refer to this array as the "stats" dictionary, which no longer contains the data

    return data, stats, hashd

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
        if data[ix2]['hash'] != data_ref[ix1]['hash']:
            pass # This is a missing particle
        else: # This is not a missing particle. It (or its duplicate) have been there.
            goodcdt[ix1] = True # Mark as  a good trajectory
            ix2 += 1 # Move onto the next data particle only if the particle wasn't missing. Otherwise, stay on this data particle.
        ix1 += 1 # Every iteration, move on to the next reference particle

    # Copy into the new array, dat3s. Could also just do a rolling update.
    data_new = np.copy(data_ref) # Make a deep copy of the reference array, as your new output array
    data_new[goodcdt] = data # Fill in all particles that were present in your data (the missing particles will stay what they were in the reference data)

    print "Fraction missing:", 1 - float(len(goodcdt[goodcdt]))/float(len(goodcdt))

    return data_new, goodcdt  # Return the data array, with all particles now present (gaps are filled)

def mpiTraj(p4dir, h5fn = None, skip=1, start=0, stop=None):
    """ Assume we have greater than one processor. Rank 0 will do the hdf5 stuff"""
    # Set some basic MPI variables
    nprocs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()
    comm = MPI.COMM_WORLD
    print "Number of processors detected through MPI:", nprocs
    
    # Check that we have at least two processors in our MPI setup    
    if nprocs < 2:
        raise Exception("Not enough processors. Need at least two, or call some other function.")

    # Everyone, get your bearings on the task to be performed
    # Get a sorted list of filenames for the pmovie files
    fns = ls.getp4(p4dir, prefix = "pmovie")[start:stop:skip] # Get list of pmovieXX.p4 files, and get the list sorted in ascending time step order
    nframes = len(fns) # Number of pmovie frames to prepare for

    # Rank0: Read in the first file and get reference data, spread that around  
    if rank == 0: # Rank 0, start working on that HDF5
        # Read the first frame, to get info like the length
        print "Rank 0 speaking, I'm reading in the first frame. Everyone else sit tight."
        data_ref, _, hashd = sortOne(fns[0])
        print "Ok, I'm going to spread that around, now."
    else:
        data_ref = None
        hashd = None # These "None" declarations are essential if we are to broadcast
    data_ref = comm.bcast(data_ref, root=0)
    hashd = comm.bcast(hashd, root=0)
    
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
        goodkeys = ['xi', 'zi', 'x', 'z', 'ux', 'uy', 'uz','q']
        # Open the HDF5 file, and step over the p4 files
        with h5py.File(h5fn, "w") as f:

            # Allocate the HDF5 datasets
            f.create_dataset("t", (nframes,), dtype='f')
            f.create_dataset("step", (nframes,), dtype='int32')
            f.create_dataset("gone", (nframes, nparts,), dtype='bool', chunks=True) # Chunking makes later retrieval along _both_ dimensions reasonably quick
            for k in goodkeys:
                f.create_dataset(k, (nframes, nparts,), dtype='f', chunks=True)
        
            # Now, iterate over the files (collect their data and save to the HDF5)
            for i in range(nframes):
                print "Rank 0: Collecting file ", i, " of ", nframes   
                t2 = dt.now()
                datdict = comm.recv(source=MPI.ANY_SOURCE, tag=i)  # Retrieve the result
                t3 = dt.now()
                datnew = datdict['datnew']
                stats = datdict['stats']
                goodcdt = datdict['goodcdt']
                badcdt = np.logical_not(goodcdt) # Flip the sign of good condit
                datnew[badcdt] = data_ref[badcdt] # Fill in the missing particles
                t4 = dt.now()

                f['t'][i] = stats['t']
                f['step'][i] = stats['step']
                f['gone'][i] = badcdt # Flag the particles that were missing
                for k in goodkeys:
                    f[k][i] = datnew[k]
                        
                t5 = dt.now()
                print "Seconds on receipt, analysis, storage: ", (t3 - t2).total_seconds(), (t4 - t3).total_seconds(), (t5 - t4).total_seconds()
                data_ref = datnew # The new array becomes the reference for next iteration
            print "ELAPSED TIME (secs)", (dt.now() - t1).total_seconds()

    else: # All other processes (non-rank 0) do the opening and reading of pmovie files, passing this info back
        print "I am a servant of the ranks."
        print "Rank", rank, ": I have frames:", myframes
        for i in range(len(myframes)):
            ix = myframes[i] # The index that refers to the entire list (of all files)
            print "Rank", rank, ": working on file", ix
            fn = fns[ix]
            dattmp, stats, _ = sortOne(fn, hashd=hashd) # Read (and sort) pmovie file
            datnew, goodcdt = fillGaps(dattmp, data_ref) # fill in missing particles (according to data_ref)
            datdict = {}
            datdict['datnew'] = datnew
            datdict['stats'] = stats
            datdict['goodcdt'] = goodcdt
            comm.send(datdict, dest=0, tag=ix) # Note: comm.isend would give an EOFError, for some reason, so don't use it.
