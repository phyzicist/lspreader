# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 08:59:43 2016

General tools used by Scott Feister (SF).

@author: Scott
"""

import os

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
    
def unChunk(l):
    """Flatten the first dimension of a list. E.g. if input is l = [[1,2,],[3,4]], output is [1,2,3,4]"""
    # Copied from http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    return [item for sublist in l for item in sublist]
    
