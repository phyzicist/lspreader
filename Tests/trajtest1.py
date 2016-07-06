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
from pmov import traj2

fns = ls.getp4(r'C:\Users\Scott\Documents\temp\mar1test\hres_osc\test', prefix='pmovie')
data_ref, stats_ref, hashd = traj2.sortOne2(fns[0])

data2, stats2, _ = traj2.sortOne2(fns[1], hashd=hashd)
data2_new = traj2.fillGaps2(data2, data_ref)