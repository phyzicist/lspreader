# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 15:30:36 2015

@author: Scott
"""

# A test of the lspreader. Specifically, reading in the second

#import h5py
import lspreader2 as rd
import matplotlib.pyplot as plt
import numpy as np
import h5py

from h5stitch2D import h5fields2D

h5path = h5fields2D(r'C:\Users\Scott\Documents\temp\lsp laser test fields')