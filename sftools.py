# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 08:59:43 2016

General tools used by Scott Feister (SF).

@author: Scott
"""

import os
import numpy as np


# NumPy arrays
def findNearest(array, value):
    """
    Get the index for the nearest value to value
    Copied from: http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    Inputs:
        array: 1D array to search for value
        value: number to search for in array
    Outputs:
        nearest value in array, index of that value
    """
    idx = np.abs(array-value).argmin()
    return array[idx], idx
    
# Vector analysis
def vecMag(a, axis=-1):
    """ Gives the magnitude of one (array of) real or complex vector """
    # a: array_like: Components of the vector
    # axis: int, optional: axis that defines the vector. By default, the last axis.
    return np.sqrt(np.sum(np.multiply(a,np.conj(a)), axis=axis)) # sqrt((a_x a_x*)+ (a_y a_y*) + ...)


# File I/O tools
def getfns(folder, ext = '', prefix = ''):
    """ Get a list of full path filenames for all files in a folder and subfolders having a certain extension
    Inputs:
       folder: a file path to the folder, e.g. "C:/hello/"
       prefix (optional): string, the prefix for the files in question, e.g. "flds_". OR: Tuple of strings, any of which can match prefix, e.g. ("flds_", "scl")
       ext (optional): string, the extension of the files in question, e.g. ".p4". OR: Tuple of strings, any of which can match suffix, e.g. (".p4", ".gz")
    Outputs:
       fns: a list of full paths to files with matching prefix and extension
    Example calls:
       pathlist1 = getfns(r'C:/myfolder1')
       pathlist1 = getfns(r'C:/myfolder2', '.png')
       pathlist2 = getfns(r'C:/myfolder3', ext=('.p4','.p4.gz'), prefix='flds')
    """
    
    fns = []
    for file in os.listdir(folder): # note 'file' can be a file, directory, etc.
        if file.endswith(ext) and file.startswith(prefix):
            fns.append(os.path.join(folder, file))
    return fns
    
def subdir(folder, name):
    """ Make a subdirectory in the specified folder, if it doesn't already exist"""
    subpath = os.path.join(folder,name)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    return subpath

# CHUNKING/UNCHUNKING TOOLS
def chunkFx(mylist, n):
    """Break list into fixed, n-sized chunks. The final element of the new list will be n-sized or less"""
    # Modified from http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    chunks = []
    for i in range(0, len(mylist), n):
        chunks.append(mylist[i:i+n])
    return chunks

def chunkIt(seq, num):
    """Divide list 'seq' into 'num' (nearly) equal-sized chunks. Returns a list of lists; some empty lists if num is larger than len(seq)."""
    # Copied directly from: http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    avg = len(seq) / float(num)
    chunks = []
    last = 0.0

    while last < len(seq):
        chunks.append(seq[int(last):int(last + avg)])
        last += avg

    return chunks

def chunkFair(mylist, nchunks):
    """ Split list into near-equal size chunks, but do it in an order like a draft pick; they all get high and low indices.
    E.g. for mylist = [0,1,2,3,4,5,6], chunkFair(mylist,4)  => [[0, 4], [1, 5], [2, 6], [3]]
    """
    chunks = [None]*nchunks
    for i in range(nchunks):
        chunks[i] = []
        
    i = 0
    while i < len(mylist):
        j = 0
        while (i < len(mylist)) & (j < nchunks):
            chunks[j].append(mylist[i])
            j += 1
            i += 1
    return chunks

def unChunk(l):
    """Flatten the first dimension of a list. E.g. if input is l = [[1,2,],[3,4]], output is [1,2,3,4]"""
    # Copied from http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    return [item for sublist in l for item in sublist]
    
