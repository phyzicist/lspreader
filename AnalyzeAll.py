# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 13:52:43 2016

Scott's analysis mission control.

Walk through a high level directory (e.g.) the data directory, and do analysis.

@author: Scott
"""

import matplotlib as mpl
#mpl.use('Agg') # Goes kind of slow! Not sure why. Something wrong with my plotting routines?

from pext import pextanalysis
import freqanalysis
import os
import re
import numpy as np

import sftools as sf
import lstools as ls
import scottplots as sp

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

def analyzeAll(p4root, outroot, pextOn=True, fldsclOn = True, freqOn=False):
    namelist = next(os.walk(p4root))[1] # Looks for all files in this list
    
    for name in namelist:
        p4dir = os.path.join(p4root, name)
        shortname = stripJunk(name) # strip the date or hyphen off the end of the folder name
        outdir = sf.subdir(outroot, shortname)
        try:
            print "Analyzing directory: " + shortname
            if pextOn:
                pextanalysis.pextFull(p4dir, outdir=outdir, shortname=shortname) # Particle extraction analysis
            if fldsclOn:
                data = ls.readFldScl(p4dir) # Read all the files involved in matching field/scalar analysis into a single data array.
                print 'Data keys:', data.keys()
                #freqanalysis.freqBatch2(data, outdir=outdir, shortname=shortname, alltime=True)
                sp.plotDens(data, outdir=outdir, shortname=shortname, alltime=True)
                sp.plotEM(data, outdir=outdir, shortname=shortname, alltime=True)
                ## Split into chunks, analyze that.
                chk_fs = 15.0#fs How many femtoseconds per (large-sized) data chunk for analysis?
                freqchunks = ls.chunkData(data, chk_fs) + ls.chunkData(data, chk_fs, offset_fs = chk_fs/2.0)
                denschunks = ls.chunkData(data, chk_fs/2.0) # Same number of chunks as bigchunks, but less long in time.
                print "Iterating over density chunks."
                for chunk in denschunks:
                    sp.plotDens(chunk, outdir=outdir, shortname=shortname)
                print "Iterating over frequency chunks."
                for chunk in freqchunks:
                    sp.plotEM(chunk, outdir=outdir, shortname=shortname)
            if freqOn:
                freqanalysis.freqFull(p4dir, outdir, nbatch = 80, divsp = 1, npool = 1) # Frequency analysis
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            # Report error and proceed
            print "Error noted. Directory: " + shortname
            pass # Move right along.
    
if __name__=='__main__':
    p4root = r'C:\Users\Scott\Documents\temp\pexttest'
    outroot = r'C:\Users\Scott\Documents\temp\pexttest\TestPlots'
    analyzeAll(p4root, outroot)
