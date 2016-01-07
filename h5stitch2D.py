# -*- coding: utf-8 -*-
"""
A library for stitching together a folder full of .p4 fields files from a 2D run

Created on Wed Dec 30 15:09:37 2015

@author: Scott
"""


import h5py
import os
import lspreader2 as rd
import numpy as np

try:
    from mpi4py import MPI
except:
    print "WARNING: MPI4PY FAILED TO LOAD. DO NOT CALL PARALLEL FUNCTIONS."
    
#except:
#    print "WARNING: Error importing mpi4py. Parallel HDF5 functions will fail."

def getfnsp4(folder, prefix = 'flds'):
    """ Get a list of full path filenames for all files 'fldsXXX.p4(.gz)' (for prefix = 'flds') in a folder, sorted by number XXX"""
    # Inputs:
    #   folder: a file path to the folder, e.g. "C:/run1", which contains files like 'fldXXX.p4'
    # Outputs:
    #   fns: a list of (sorted) full paths to files matching 'fldsXXX.p4', e.g. fns = ['C:/run1/flds1.p4', 'C:/run1/flds5.p4', 'C:/run1/flds10.p4']
    
    fns = [] # Will store the list of filenames
    nums = [] # Will store a list of the XXX numbers in 'fldsXXX.p4'
    
    pattern = r'^(' + prefix + ')(\d*?)(\.p4)(\.gz){0,1}$'
    # If prefix = 'flds', matches 'fldsXXX.p4' and 'fldsXXX.p4.gz', where XXX is a number of any length. Does not match 'dfldsXXX.p4' or 'fldsXXX.p4.zip'
    

    for name in os.listdir(folder): # note 'file' can be a file, directory, etc.
        m = re.search(pattern, name)
        if m: # If the filename matches the pattern, add this file to the list
            fns.append(os.path.join(folder, name))
            
            # Extract "XXX" as a number from "fldsXXX.p4(.gz)"
            # Note: re.split(r'^(fld)(\d*?)(\.p4)(\.gz){0,1}$', 'flds5020.p4.gz) returns something like ['', 'fld', '5020', '.p4', '.gz', ''], so we take the element at index 2 and convert to integer
            fnum = int(re.split(pattern, name)[2])
            nums.append(fnum)

    # Sort the field names by their XXX number    
    idx = np.argsort(nums) # Get a list of sorted indices such that the filenames will be sorted correctly by time step
    fns = np.array(fns)[idx].tolist() # Convert filenames list to a numpy string array for advanced indexing, then convert back to a python list
    return fns

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

def h5fields2Dser(folder, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], flds = ['E', 'B']):
    """ 100% serial processing. Read in a folder full of flds*.p4 file, stitching them together assuming 2D assumptions, and create an HDF5 file """
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