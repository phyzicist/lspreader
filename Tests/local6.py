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
import matplotlib.pyplot as plt
import numpy as np
import sftools as sf
from pext import pextanalysis as pa
import scipy.constants as sc
from pext import quantities
    
#p4root = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump'
#outroot = r'C:\Users\Scott\Documents\temp\sclrtest\analysis'
#aa.analyzeSome(p4root, outroot, ['a50f-14_mres_so'])


#p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test0-2015-15-21_1923'
#p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test2'
#data = ls.readFldScl(p4dir)
#
#shortname = 'test0'
#outdir = r'C:\Users\Scott\Documents\temp\sclrtest\analysis\test0'
#pa.pextFull(p4dir, outdir=outdir, shortname=shortname)
#sp.plotBScat(data, tcut=65, outdir=outdir, shortname = shortname)


#chk_fs = 1.0#fs How many femtoseconds per (large-sized) data chunk for analysis?
#freqchunks = ls.chunkData(data, chk_fs) + ls.chunkData(data, chk_fs, offset_fs = chk_fs/2.0)
#denschunks = ls.chunkData(data, chk_fs/2.0) # Same number of chunks as bigchunks, but less long in time.


### DENSITY VISUAL TEST
p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\byc1f-6_mres_so'
data = ls.readFldScl(p4dir)
sp.plotDens(data, outdir=p4dir, shortname = '', alltime=False)

xgv = data['xgv']*1e4 # x values in microns
zgv = data['zgv']*1e4
dx = np.mean(np.diff(xgv))# dx in microns
dz = np.mean(np.diff(zgv))

### CALCULATIONS
## Mean electron density
#edens = np.mean(data['RhoN10'],0)

## PEXTTEXT
#p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\a28f-14_mres_so'
#outdir = r'C:\Users\Scott\Documents\temp\sclrtest\analysis\a28f-14_mres_soTEST'
#pa.pextFull(p4dir, outdir = outdir, shortname = 'a28f-14_mres_soTEST', Utot_Jcm = 25)
