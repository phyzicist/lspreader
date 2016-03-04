# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 17:21:46 2016

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

def myhash(dat):
    #h = dat['xi'].astype('f8') * (1e2 + dat['yi'].astype('f8')) * (1e4*dat['zi'].astype('f8'))
    h = dat['xi']
    #h = dat['xi'] + dat['yi'] + dat['zi']
    return h

dat1s = np.sort(dat, order=['xi','yi','zi'])

allhash = myhash(dat1s)

## FIRST ROUND: Separate the good and bad hashes
_, goodix = np.unique(allhash, return_index=True)
badix = np.setdiff1d(range(len(allhash)), goodix)
goodhash = allhash[goodix]
badhash = np.unique(allhash[badix])

## SECOND ROUND : order by hash, fill in the gaps
#STRIP OUT THE BAD HASHES FROM NUMPY ARRAY
dat2s = np.sort(dat, order=['xi','yi','zi'])

allhash2 = myhash(dat2s)
goodcdt = np.zeros(allhash2.shape, dtype=bool)

ix_hash = 0
goodhashp = np.concatenate((goodhash, [0]))
for i in range(len(allhash2)):
    if allhash2[i] == goodhashp[ix_hash]:
        goodcdt[i] = True
        ix_hash += 1

#_, goodix = np.setdiff1d(allhash, badhash, return_index=True)
#dat1s[allhash!=badhash]

print "Non-unique fraction:", len(nonunique)/float(len(unique_ix))

## PMOVIE TEST
fn = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc\pmovie1.p4.gz"
frames = rd.read_movie2(fn)
dat = frames[0]['data']

myhashes = myhash(dat)

# FIRST DATASET (SIM)
subdat1 = dat[:5]

print np.round(990-subdat1['xi']*1e6)

subdat1s = np.sort(subdat1, order=['xi','yi','zi'])
hash1 = myhash(subdat1)

print np.round(990-subdat1s['xi']*1e6)

# SECOND DATSET (SIM)
subdat2 = dat[[1,2,4]]
subdat2s = np.sort(subdat2, order=['xi','yi','zi'])
hash2 = subdat2['xi']

print np.round(990-subdat2s['xi']*1e6)

print np.setdiff1d(hash1, hash2)