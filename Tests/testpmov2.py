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
    #h = dat['xi']
    #h = dat['xi'] + dat['yi'] + dat['zi']
    h = dat['xi']
    return h


## PMOVIE TEST
fn = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc\pmovie1.p4.gz"
frames = rd.read_movie2(fn)
dat = frames[0]['data']

dat1s = np.sort(dat, order=['xi','yi','zi'])

allhash = myhash(dat1s)

## FIRST ROUND: Separate the good and bad hashes
_, goodix = np.unique(allhash, return_index=True)
badix = np.setdiff1d(range(len(allhash)), goodix)
goodhash = allhash[goodix]
badhash = np.unique(allhash[badix])

dathash = dat1s[goodix]


## SECOND ROUND : order by hash, fill in the gaps
#STRIP OUT THE BAD HASHES FROM NUMPY ARRAY
dat2s = np.sort(dat, order=['xi','yi','zi'])

allhash2 = myhash(dat2s)
goodcdt = np.zeros(allhash2.shape, dtype=bool)

ix_hash = 0
#goodhashp = np.concatenate((goodhash, [0]))
#for i in range(len(allhash2)):
#    if allhash2[i] == goodhashp[ix_hash]:
#        goodcdt[i] = True
#        ix_hash += 1

goodcdt = np.zeros(dat1s.shape, dtype=bool)
ix1 = 0
ix2 = 0
while ix2 < len(dat2s):
    # Check if any xi, zi, or yi are unequal; if not, label this as a "good one".
    if (dat2s[ix2]['xi'] != dat1s[ix1]['xi']) & (dat2s[ix2]['zi'] != dat1s[ix1]['zi']) & (dat2s[ix2]['yi'] != dat1s[ix1]['yi']):
        pass
    else:
        goodcdt[ix1] = True
        ix2 += 1
    ix1 += 1

#_, goodix = np.setdiff1d(allhash, badhash, return_index=True)
#dat1s[allhash!=badhash]

print "Non-unique fraction:", len(goodix)/float(len(badix))
print "Number of unique hashes:", len(goodix)
print "Number of matches (should be equal):", len(goodcdt[goodcdt])
print "MISSED MATCHES:", len(goodix) - len(goodcdt[goodcdt])