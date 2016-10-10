#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
chargedist_cylinder.py: Test problem exploring how quickly I can calculate electric potential of a cylinder of charge

Created by Scott Feister on Mon Oct 10 16:23:43 2016

Will take a 2D profile of a cylinder, but know that it extends in the Z dimension.
Will start with cylinder of large but finite length 2a, and radius R.

Reference for solution: https://www.miniphysics.com/uy1-electric-potential-of-an-infinite-line-charge.html
"""

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print "Hello."
    
    a = 1.0e4 # Set a as one cm, 1e4 microns
    rc = 3.0 # Set cylinder radius as 3 microns
    
    # Create the cylinder charge background
    #l = # Line charge density (coulombs/micron) ### CONTINUE WORK HERE. -SF 10-10-2016
    nx = ny = 100 # Number of points in each dimension
    xgv = np.linspace(-30, 30, nx)
    ygv = np.linspace(-30, 30, ny)
    X, Y = np.meshgrid(xgv, ygv)
    Rcent = np.sqrt(X**2 + Y**2)
    pass
