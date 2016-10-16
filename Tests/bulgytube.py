#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
name.py: description

Created by Scott Feister on Fri Oct 14 21:42:19 2016
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import math
import numpy as np

    
if __name__ == "__main__":    
    zmin=-10
    zmax=10
    xmin=-20
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    D = np.zeros_like(X)
    
    mu = -5
    sigma = 3
    ct = (np.abs(Z) < 1 + 10*mlab.normpdf(X, mu, sigma))
    D[ct] = 1.0
    clabel="Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)
    