# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 15:30:36 2015

@author: Scott
"""

# A test of the lspreader. Specifically, reading in the second

#import h5py
import lspreader2 as rd
import matplotlib.pyplot as plt
import numpy as np
import h5py

import os

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
    
fn = r"C:\Users\Scott\Documents\temp\lsp laser test fields\flds2090.p4"

fns = getfns(r'C:\Users\Scott\Documents\temp\lsp laser test fields',ext = '.p4', prefix = 'flds')

#for i in range(len(fns)):
#    h5f = h5py.File('data.h5', 'w')
#    h5f.create_dataset(str(i), data=doms)

# Read in this timestep
doms, header = rd.read_flds2(fn, flds=['E','B'])
# Extract the interesting fields and stitch together for a template
nfiles = len(fns)
fld_ids = ['Ex','Ey','Ez','Bx','By','Bz']

# Open the HDF5 file
with h5py.File('data3.h5','w') as f:
    # Pre-allocate the HDF5 arrays
    _, xgv, zgv = rd.stitch2D(doms, fld_ids[0])
    
    f.create_dataset('times', (nfiles,), compression='gzip', dtype='float32')
    f.create_dataset('xgv', compression='gzip', data = xgv)
    f.create_dataset('zgv', compression='gzip', data = zgv)
    f.create_dataset('filenames', data = fns) # A bit overkill to store the full filenames, but who cares? Space is not too bad.
    for k in fld_ids:
        f.create_dataset(k,(nfiles,len(zgv),len(xgv)), compression='gzip', dtype='float32')
    
    # Read in the data and fill up the datafile
    for i in range(nfiles): # Iterate over the files
        fn = fns[i]
        doms, header = rd.read_flds2(fn, flds=['E','B'])
        f['times'][i] = header['timestamp']
        for fld_id in fld_ids: # Iterate over the requested fields, stitching then adding them to HDF5 file
            fld, xgv, zgv = rd.stitch2D(doms, fld_id)
            f[fld_id][i,:,:] = fld
    
#
#for i in range(len(fld_ids)):
#    fld_id = fld_ids[i]
#    fld, xgv, zgv = rd.stitch2D(doms, fld_id)
#
#
#compress=False
#f = h5py.File('data.hdf5','w')
#for i in range(len(doms)):
#    grp = f.create_group(str(i))
#    d = doms[i]
#    for k in d:
#        if k not in grp.keys():
#            if compress:
#                grp.create_dataset(k,data=d[k],compression="lzf")
#            else:
#                grp.create_dataset(k,data=d[k])
#        else:
#            grp[k] = d[k]
#f.close()