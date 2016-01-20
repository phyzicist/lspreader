# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:57:20 2016

@author: Scott
"""

import sys
try:
    import lspreader2 as rd
except:
    print "Modifying path to include LSPreader"
    readerpath = r'C:\Users\Scott\Documents\Programming\Python\lspreader'
    sys.path.append(readerpath) # Add the LSP reader as taking precedence after all but the current directory
    import lspreader2 as rd

import scottplots as sp
import lstools as ls
import AnalyzeAll as aa

p4root = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump'
outroot = r'C:\Users\Scott\Documents\temp\sclrtest\analysis'
aa.analyzeAll(p4root, outroot)


#p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test0'
#data = ls.readFldScl(p4dir)
#
#shortname = 'test0'
#outdir = r'C:\Users\Scott\Documents\temp\sclrtest\analysis\test0'
#sp.plotme(data, outdir= outdir, shortname=shortname)

#chk_fs = 1.0#fs How many femtoseconds per (large-sized) data chunk for analysis?
#freqchunks = ls.chunkData(data, chk_fs) + ls.chunkData(data, chk_fs, offset_fs = chk_fs/2.0)
#denschunks = ls.chunkData(data, chk_fs/2.0) # Same number of chunks as bigchunks, but less long in time.
