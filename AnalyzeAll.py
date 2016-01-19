# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 13:52:43 2016

Scott's analysis mission control.

Walk through a high level directory (e.g.) the data directory, and do analysis.

@author: Scott
"""

import matplotlib as mpl
#mpl.use('Agg') # Goes kind of slow! Not sure why. Something wrong with my plotting routines?
import sys

try:
    import lspreader2 as rd
except:
    print "Modifying path to include LSPreader"
    readerpath = r'C:\Users\Scott\Documents\Programming\Python\lspreader'
    sys.path.append(readerpath) # Add the LSP reader as taking precedence after all but the current directory
    import lspreader2 as rd

from pext import pextanalysis
import os
import re

def subdir(folder, name):
    """ Make a subdirectory in the specified folder, if it doesn't already exist"""
    subpath = os.path.join(folder,name)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    return subpath

def stripJunk(name):
    """ Strips a date string off the end of a folder name, giving a "shorter" version of the name.
    Very specific to the way I've been saving directories as dirname = shortname + date. I don't want the date, though.
    """
    m1 = re.match(r'(.*?)(-\d{4}-\d{2}-\d{2})(.*)', name) # e.g. 'a50f-14xL_mres_so-2016-01-17_2006' => 'a50f-14xL_mres_so'
    m2 = re.match(r'(.*?)(-)$', name) # e.g. '2d-nosolid-' => '2d-nosolid'
    if m1:
        shortname = m1.groups()[0]
    elif m2:
        shortname = m2.groups()[0]
    else:
        shortname = name
    return shortname

if __name__=='__main__':
    p4root = r'C:\Users\Scott\Documents\temp\pexttest'
    outroot = r'C:\Users\Scott\Documents\temp\pexttest\TestPlots'
    namelist = next(os.walk(p4root))[1] # Looks for all files in this list
    
    for name in namelist:
        p4dir = os.path.join(p4root, name)
        shortname = stripJunk(name)
        outdir = subdir(outroot, shortname)
        try:
            print "Analyzing directory: " + shortname
            pextanalysis.fullAnalyze(p4dir, outdir=outdir, shortname=shortname)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            # Report error and proceed
            print "Error noted. Directory: " + shortname
            pass # Move right along.
