# -*- coding: utf-8 -*-
"""
Created on Thu Jan 07 11:30:19 2016

@author: Scott
Sorting the files

"""
# Documentation for reg exp of strings at: https://docs.python.org/2/library/re.html

# A regular expression (or RE) specifies a set of strings that matches it; the functions in this module let you check if a particular string matches a given regular expression


import re
import os
import numpy as np

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

fns = getfnsp4(r'C:\Users\Scott\Documents\temp\lastest\run2_lastest longer\gzipdir', prefix = 'flds')
print fns[0]