# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 18:36:07 2016

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
import matplotlib as mpl
from pmov import traj2

if __name__ == '__main__':
    # TODO: Save the dat3s slice into the HDF5 file
    mydir = r"C:\Users\Scott\Documents\temp\mar1test\hres_osc"
    #dat1 = traj2.parTraj(mydir, nprocs=8)
    traj2.mpiTraj(mydir)
