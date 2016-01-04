# -*- coding: utf-8 -*-
"""
A library for stitching together a folder full of .p4 fields files from a 2D run

Created on Wed Dec 30 15:09:37 2015

@author: Scott
"""


import h5py
import os
import lspreader2 as rd

def getfns(folder, ext = '', prefix = ''):
    """ Get a list of full path filenames for all files in a folder and subfolders having a certain extension"""
    # Inputs:
    #   folder: a file path to the folder, e.g. "C:/hello/"
    #   prefix (optional): the prefix for the files in question, e.g. "flds_"
    #   ext (optional): the extension of the files in question, e.g. ".p4"
    # Outputs:
    #   fns: a list of full paths to files with matching prefix and extension
    
    fns = []
    for file in os.listdir(folder): # note 'file' can be a file, directory, etc.
        if file.endswith(ext) and file.startswith(prefix):
            fns.append(os.path.join(folder, file))
    return fns



def chunkIt(seq, num):
    """Divide list 'seq' into 'num' (nearly) equal-sized chunks. Returns a list of lists; some empty lists if num is larger than len(seq)."""
    # Copied directly from: http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def h5fields2D(folder, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], flds = ['E', 'B']):
    """ Read in a folder full of flds*.p4 file, stitching them together assuming 2D assumptions, and create an HDF5 file """
    if not h5path:
        h5path = os.path.join(folder, 'fields2D.hdf5')

    fns = getfns(folder, ext = '.p4', prefix = 'flds')
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

def h5fields2Da(folder, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], flds = ['E', 'B']):
    """ Step A of: Read in a folder full of flds*.p4 file, stitching them together assuming 2D assumptions, and create an HDF5 file """
    # In Step A (to be run on a single processor before any other activity), we make and allocate the HDF5 file.
    if not h5path:
        h5path = os.path.join(folder, 'fields2D.hdf5')

    fns = getfns(folder, ext = '.p4', prefix = 'flds')
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

    return h5path

def h5fields2Db(folder, nprocs = 1, rank = 0, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], flds = ['E', 'B']):
    """ Step B of: Read in a folder full of flds*.p4 file, stitching them together assuming 2D assumptions, and create an HDF5 file """
    # In Step B (which can be called by many processors, or just one by leaving the defaults), we append the data to the newly-created HDF5 file.

    if not h5path:
        h5path = os.path.join(folder, 'fields2D.hdf5')

    fns = getfns(folder, ext = '.p4', prefix = 'flds') # the COMPLETE list of filenames (for all processors)
    nfiles = len(fns) # the total number of files (for all processors)

    idx_all = range(nfiles) # This list of indices can be split up later across processors, outside this function
    idx_chunked = chunkIt(idx_all, nprocs) # Break the list into N chunks, where N = number of processors
    idx_part = idx_chunked[rank] # Assign each processor one chunk of the list. Therefore, idx_part (list of indices) specifies which elements to analyze here.

    with h5py.File(h5path,'r+') as f:
        # Read in all the data and fill up the HDF5 file
        for i in idx_part: # Iterate over the indices given for this processor
            fn = fns[i]
            doms, header = rd.read_flds2(fn, flds=flds)
            f['times'][i] = header['timestamp']
            for fld_id in fld_ids: # Iterate over the requested fields, stitching then adding them to HDF5 file
                fld, xgv, zgv = rd.stitch2D(doms, fld_id)
                f[fld_id][i,:,:] = fld
    return h5path, idx_part
