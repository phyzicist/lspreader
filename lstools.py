# -*- coding: utf-8 -*-
"""
A library with somewhat general tools for 2D LSP analysis of the type run for the electron experiment. -Scott

Renamed from h5stitch2D.py.

Created on Wed Dec 30 15:09:37 2015

Changelog:
    2015-01-09: Updated the way the 'x' is stripped off 'Ex' to work with scalars, too
    
@author: Scott
"""


import h5py
import os
import lspreader2 as rd
import numpy as np
import re
import gzip
import re

try:
    from mpi4py import MPI
except:
    print "WARNING: MPI4PY FAILED TO LOAD. DO NOT CALL PARALLEL FUNCTIONS."

def getp4(folder, prefix = 'flds'):
    """ Get a sorted list of full path filenames for all files 'fldsXXX.p4(.gz)' (for prefix = 'flds') in a folder, sorted by number XXX"""
    # Inputs:
    #   folder: a file path to the folder, e.g. "C:/run1", which contains files like 'fldXXX.p4'
    #   prefix: string with which the filename starts, e.g. 'flds' for files matching '.../flds*.p4'
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


def getFldScl(p4dir):
    """ Get a 1-to-1 list of fields and scalars filenames in a directory; that is, exclude any fld or scl files that are missing their partner. """
    fld = getp4(p4dir, 'flds')
    scl = getp4(p4dir, 'sclr')
    
    fns_fld = []
    fns_scl = []

    for fn_scl in scl:
        folder, name = os.path.split(fn_scl)
        name_fld = name.replace('sclr','flds')
        fn_fld = os.path.join(folder, name_fld)
        if fn_fld in fld:
            fns_fld.append(fn_fld)
            fns_scl.append(fn_scl)
    
    return fns_fld, fns_scl

def readFldScl(p4dir, divsp=1):
    """ Read matching field and scalar files (default fld_ids) in a directory into a single data array """
    ## READ MATCHING FIELDS AND SCALARS INTO DATA ARRAY
    fns_fld, fns_scl = getFldScl(p4dir)
    
    data_fld = fields2D(fns_fld, divsp=divsp)
    data_scl = scalars2D(fns_scl, divsp=divsp)
    
    # Sanity check that times, xgv, and zgv are essentially identical between the two
    if np.max(np.abs(data_scl['xgv'] - data_fld['xgv'])) + np.max(np.abs(data_scl['zgv'] - data_fld['zgv'])) + np.max(np.abs(data_scl['times'] - data_fld['times'])) > 0.000001:
        raise Exception("Scalar and field files grid vectors (or times) do not match!! This shouldn't be the case. Major problem with file choice or with lspreader.")
    
    data = {} # The merged data array. Note that no copies of arrays will be made, just pointers shuffled around.
    for k in data_fld.keys():
        data[k] = data_fld[k]
    for k in set(data_scl.keys()).difference(data_fld.keys()): # Add only those keys not already in the key list
        data[k] = data_scl[k]
    
    return data
    
def subData(data, startwith, stopbefore, skip = 1):
    """ Get a subset of the data dict, along the time axis. (Not a deep copy; just shuffling pointers. That is: you change data_sub elements, you change data.)
    Inputs:
        data: dictionary of numpy arrays, of the form of the output from fields2D()
        startwith, stopbefore, skip: The indices with which to slice time. If 'times' were a NumPy array:
            then what you get out is eqivalent to times[startwith:stopbefore:skip]. By default, skip over no data (skip=1).
    Outputs:
        data_sub: array of the same form as data, but smaller along the time axis. Now has time length 'stop - start'
    Example call:
        data = fields2D(filenames) # data dict is as many times steps long as there are filenames
        data_sub1 = subData(data, 0, 3, skip = 1) # data_sub1 dict contains only steps 0, 1, and 2
        data_sub2 = subData(data, 1, 11, skip = 3) # data_sub2 dict contains only steps 1, 4, 7, and 10
    """    
    data_sub = {}
    for k in data.keys():
        if k in set(['xgv','ygv','zgv']):
            data_sub[k] = data[k] # No time axis exists in these arrays. Place into data_sub without modification.
        else:
            data_sub[k] = data[k][startwith:stopbefore:skip] # Get a subset along the time axis (axis 0), then place into data_sub.
    return data_sub

def chunkData(data, chk_fs, offset_fs = 0):
    """Break data along the time step axis into fixed, n-sized chunks (where n is determined by desired chunk size in fs).
    The final element of the new list will be n-sized as well; some leftover frames at the end will not be included
    
    Inputs:
        data: the usual data dict (which is N frames long)
        chk_fs: desired chunk size, in femtoseconds
        offset_fs: skip this many fs past the 0th frame before starting the chunking
    Outputs:
        chunks: list of data dicts  (each being n frames long, where n <= N)       
    """
    # Modified from http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    dt_fs = np.mean(np.diff(data['times']*1e6)) # Femtoseconds between frames in data (assumed fixed)
    n = int(np.ceil(chk_fs / dt_fs)) # Desired # frames per chunk, based on desired fs per chunk
    N = len(data['times']) # Total number frames in the data dict
    offset = max(0, int(np.round(offset_fs / dt_fs))) # by default, offset is 0, meaning we start chunking from the 0th frame.
    chunks = [] # The output will be a list of data dicts
    for i in range(offset, N + 1 - n, n):
        chunk = subData(data, i, i + n)
        chunks.append(chunk)
    return chunks

def fields2D(fns, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], divsp=1, pool = None):
    """ Read in the flds*.p4(.gz) files in the list fns, stitching them together assuming 2D assumptions, and create output arrays. The read-in occurs in parallel if the input parameter pool is set (Pool of multiprocessing threads e.g. via Pool(10)). 
    Inputs:        
        divsp: integer, divisor by which to reduce the spatial resolution (e.g. divsp = 2 reduces field dimensions from 300x200 to 150x100)
    Changelog:
        2016-01-19 Should also work with scalar input filenames, as needed.
    """
    # The first dimension of each array must be nfiles long.
    
    ## Extract "E" from "Ex" in fld_ids (for rd.read_flds2() call later)    
    flds = list(set([re.sub('[xyz]$','',s) for s in fld_ids])) # we strip the "x","y","z" last character off our fields, then set() gives only unique elements of a list, and list() converts this set back to list
    # E.g. if fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], flds = ['E','B']
    # Alternatively, if fld_ids = ['RhoxN1', 'RhoxN2'], flds = ['RhoxN1', 'RhoxN2']
        
    ## Read in the first file as a template for data array pre-allocations
    doms, header = rd.read_flds2(fns[0], flds=flds)
        
    nfiles = len(fns) # Count the number of files we need to read
    
    ## Pre-allocate the NumPy arrays inside an output dict called 'data'
    data = {} # define 'data' as a python dictionary that will store all the data read in, including fields, as NumPy arrays    
    _, xgv, zgv = rd.stitch2D(doms, fld_ids[0], divsp = divsp) # Extract the interesting fields and stitch together for a template
    data['times'] = np.zeros((nfiles,))
    data['xgv'] = xgv
    data['zgv'] = zgv
    data['filenames'] = np.array(fns)

    for k in fld_ids:
        data[k] = np.zeros((nfiles,len(zgv),len(xgv))) # The vector field elements

    ## Read in the files
    if pool: # OPTION A: PARALLEL READ OF FILES INTO DATA DICT
        print "Using parallel pool to read", len(fns), "files. (This could take a little while.)"
        args_iter = zip(fns, [fld_ids]*nfiles, [flds]*nfiles, [divsp]*nfiles) # Make a list nfiles long, with tuples of (filename, fld_ids, flds)
        outs = pool.map(readOne, args_iter)
        flds1, times = zip(*outs)
        flds1 = np.array(flds1)
        data['times'] = np.array(times) # Save the timestamps into the data dict
        # Re-order and save into data
        for i in range(len(fld_ids)):
            k = fld_ids[i]
            data[k] = flds1[:,i,:,:] # Save the field elements into the data dict
        print "Done reading files."

    else: # OPTION B: SERIAL READ OF FILES INTO DATA DICT
        print "Reading the files in serial."        
        # Read in all the data and fill up the NumPy arrays
        for i in range(nfiles): # Iterate over the files
            print "Reading file", i, "of", nfiles        
            fn = fns[i]
            doms, header = rd.read_flds2(fn, flds=flds)
            data['times'][i] = header['timestamp']
            for k in fld_ids: # Iterate over the requested fields, stitching then adding them to fldDict arrays
                fld, _, _ = rd.stitch2D(doms, k, divsp = divsp)
                data[k][i,:,:] = fld

    return data

def scalars2D(fns, fld_ids = ['RhoN1', 'RhoN2', 'RhoN3',  'RhoN4', 'RhoN5', 'RhoN6', 'RhoN7', 'RhoN8', 'RhoN9', 'RhoN10',  'RhoN11'], divsp=1, pool = None): 
    """ A wrapper for fields2D, with default inputs better suited for scalar filenames input."""
    # ['RhoN1', 'RhoN2', 'RhoN3',  'RhoN4', 'RhoN5', 'RhoN6', 'RhoN7', 'RhoN8', 'RhoN9', 'RhoN10',  'RhoN11'] Each corresponding to density of species 1, 2, 3..., 11
    # ['RhoN2', 'RhoN10',  'RhoN11']  # oxygen+, electrons, protons in our sims  
    return fields2D(fns, fld_ids = fld_ids, divsp=1, pool = pool)
    
def readOne(args):
    """ Helper function for multiprocessing Pool.map() call of Parallel read for fields2D. Reads in the fields from one file."""
    # Pool.map cannot have more than one input; hence, we have an input that is tuple of (filename, fld_ids, flds).
    # Output is a tuple of read-in fields (for all fields in fld_ids) and timestamp.
    fn = args[0]
    fld_ids = args[1]
    flds_label = args[2]
    divsp = args[3]

    doms, header = rd.read_flds2(fn, flds=flds_label)

    time = header['timestamp']
    for i in range(len(fld_ids)): # Iterate over the requested fields, stitching then adding them to fldDict arrays
        fld, _, _ = rd.stitch2D(doms, fld_ids[i], divsp = divsp)
        if i == 0:
            flds1 = np.zeros((len(fld_ids), fld.shape[0], fld.shape[1]))
            flds1[i,:,:] = fld
        else:
            flds1[i,:,:] = fld
    return (flds1, time)


def getTimes(fns):
    """ Get a list of times (in ns), from a list of .p4(.gz) filenames. """

    nfiles = len(fns)
    times = np.zeros(nfiles)    
    for i in range(nfiles):
        fn = fns[i]
        _, ext = os.path.splitext(fn)
        
        if ext == '.gz':
            with gzip.open(fn, 'rb') as f:
                header = rd.get_header(f)
        else:
            with open(fn, 'rb') as f:
                header = rd.get_header(f)
        times[i] = header['timestamp']

    return times

def h5fields2D(folder, h5path=None, fld_ids = ['Ex','Ey','Ez','Bx','By','Bz'], pool = None):
    """ Serial or parallel read-in of files. Works great. Creates the HDF5 file without all the pain. """
    if not h5path:
        h5path = os.path.join(folder, 'fields2D.hdf5')
    fns = getp4(folder)
    nfiles = len(fns)
    print 'Total number of files:', len(fns)

    print "Opening the HDF5 file"
    with h5py.File(h5path,'w') as f:
        # Read all the files into RAM
        print "Reading files into NumPy arrays in RAM"
        data = fields2D(fns, fld_ids=fld_ids, pool=pool)
        # Build the HDF5 file, assuming every element in "data" is a NumPy array
        print "Saving arrays in RAM to HDF5"
        for k in data:
            f.create_dataset(k, data = data[k], compression='gzip', compression_opts=4)
    print "All done!"
    return h5path
