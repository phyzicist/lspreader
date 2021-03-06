#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
loontunes.py: Wild exploration of target designs

Created by Scott Feister on Fri Oct 14 20:17:18 2016
"""

from threescale import write2DLSP
from plaswriter import critDens

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.ndimage.filters import gaussian_filter

def ccutwire():
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
    
    x_crit = -15
    x_right = -14
    x_right2 = -8
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_crit + 1) & (np.abs(Z) < 3.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[D < 0.5*n_crit] = 0.0
    
    scale2 = 15
#    ct = (X >= x_crit + 1) & (X < x_right2) & (np.abs(Z) < 0.4*(X - 2*x_crit))
    ct = (X >= x_crit + 1) & (X < x_right2) & (np.abs(Z) < 0.9*(-X - 7))
    #ct = (X >= x_right) & (X < x_right2) & (np.abs(Z) < 1.)
    D[ct] = 0.3*n_crit * np.exp(-((X[ct] - x_crit))/scale2)

    scale2 = 15
    ct = (X >= x_right2) & (np.abs(Z) < 1.)
    D[ct] = 0.3*n_crit * np.exp(-((X[ct] - x_crit))/scale2)


    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\ccutwire.dat"
    write2DLSP(fname, xgv, zgv, D)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def streamwire():
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
    
    x_crit = -15
    x_right = -14
    x_right2 = -8
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_crit + 1) & (np.abs(Z) < 3.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[D < 0.5*n_crit] = 0.0
    
    scale2 = 7.2
    mu = -8
    sigma = 4
    ct = (X >= x_right) & (np.abs(Z) < 1 + 10*mlab.normpdf(X, mu, sigma))
    D[ct] = 0.5*n_crit * np.exp(-((X[ct] - x_crit))/scale2)
    

    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\streamwire.dat"
    write2DLSP(fname, xgv, zgv, D)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def gridwire():
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
    
    x_crit = -15
    x_right = -14
    x_right2 = -8
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_crit + 1) & (np.abs(Z) < 7.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    
    scale2 = 15
#    ct = (X >= x_crit + 1) & (X < x_right2) & (np.abs(Z) < 0.4*(X - 2*x_crit))
    ct = (X >= x_crit + 1) & (X < x_right2) & (np.abs(Z) < 0.9*(-X - 7))
    #ct = (X >= x_right) & (X < x_right2) & (np.abs(Z) < 1.)
    D[ct] = 0.2*n_crit * np.exp(-((X[ct] - x_crit))/scale2)

    scale2 = 15
    ct = ((X >= x_right2) & (np.abs(Z) < 1.)) | ((X > -6) & (X < -5)) | ((X > -3) & (X < -2)) | ((X > 0) & (X < 1)) | ((X > 3) & (X < 4)) | ((X > 6) & (X < 7)) | ((X > 8.7) & (X < 9.7)) | ((X > -8) & (np.abs(Z) > 4.0))
    D[ct] = 0.2*n_crit * np.exp(-((X[ct] - x_crit))/scale2)


    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\gridwire.dat"
    write2DLSP(fname, xgv, zgv, D)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def grid2wire():
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
    
    x_crit = -15
    x_right = -14
    x_right2 = -12
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_crit + 1) & (np.abs(Z) < 7.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    
    scale2 = 15
#    ct = (X >= x_crit + 1) & (X < x_right2) & (np.abs(Z) < 0.4*(X - 2*x_crit))
    ct = (X >= x_crit + 1) & (X < x_right2) & (np.abs(Z) < 0.9*(-X - 11))
    #ct = (X >= x_right) & (X < x_right2) & (np.abs(Z) < 1.)
    D[ct] = 0.2*n_crit * np.exp(-((X[ct] - x_crit))/scale2)

    scale2 = 15
    ct = ((X >= x_right2) & (np.abs(Z) < 1.)) | ((X > -11) & (X < -10)) | ((X > -8.7) & (X < -7.6)) | ((X > -6) & (X < -5)) | ((X > -3) & (X < -2)) | ((X > 0) & (X < 1)) | ((X > 3) & (X < 4)) | ((X > 6) & (X < 7)) | ((X > 8.7) & (X < 9.7)) | ((X > -12) & (np.abs(Z) > 4.0))
    D[ct] = 0.15*n_crit


    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\grid2wire.dat"
    write2DLSP(fname, xgv, zgv, D)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def stream2wire():
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
    
    x_crit = -15
    x_right = -14
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_crit + 1) & (np.abs(Z) < 4.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[D < 0.25*n_crit] = 0.0
    
    scale2 = 30
    mu = -8
    sigma = 4
    ct = (X >= x_right) & (np.abs(Z) < 0.5 + 10*mlab.normpdf(X, mu, sigma))
    D[ct] = 0.25*n_crit *np.exp(-(-(X[ct] + 10))/scale2)
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stream2wire.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def stream3wire():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    
    x_crit = -17
    x_right = -14
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_right) & (np.abs(Z) < 4.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[D < 0.2*n_crit] = 0.0
    
    scale2 = 30
    mu = -8
    sigma = 4
    ct = (X >= x_right) & (np.abs(Z) < 0.5 + 10*mlab.normpdf(X, mu, sigma))
    D[ct] = 0.25*n_crit *np.exp(-(-(X[ct] + 10))/scale2)
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stream3wire.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def stream5wire():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    
    x_crit = -17
    x_right = -14
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_right) & (np.abs(Z) < 8.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[(D < 0.2*n_crit) & (X < x_right)] = 0.2*n_crit
    
    scale2 = 30
    mu = -8
    sigma = 4
    ct = (X >= x_right) & (np.abs(Z) < 0.5 + 10*mlab.normpdf(X, mu, sigma))
    D[ct] = 0.8*n_crit
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stream5wire.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def stream6wire():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    
    x_crit = -17
    x_right = -14
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_right) & (np.abs(Z) < 5.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[(D < 0.2*n_crit) & ct] = 0.2*n_crit
    
    scale2 = 30
    mu = -8
    sigma = 4
    ct = (X >= x_right) & (np.abs(Z) < 0.5 + 10*mlab.normpdf(X, mu, sigma))
    D[ct] = 0.8*n_crit
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stream6wire.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def stream7wire():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    
    x_crit = -17
    x_right = -14
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_right) & (np.abs(Z) < 5.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[(D < 0.2*n_crit) & ct] = 0.0
    
    scale2 = 30
    mu = -8
    sigma = 4
    ct = (X >= x_right) & (np.abs(Z) < 0.5 + 10*mlab.normpdf(X, mu, sigma))
    D[ct] = 0.8*n_crit
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stream7wire.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def slabcomb():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    
    x_crit = -17
    x_right = -14
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_right) & (np.abs(Z) < 5.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[(D < 0.2*n_crit) & ct] = 0.2*n_crit
    

    # Middle bands
    zbands = np.arange(-10, 10, 1) # Location of bands
    bandwid = 0.2 # Width of bands
    
    x_sink = 1.0 # Location of the current sink
    ct1 = np.min(np.abs(np.tile(Z, (len(zbands),1,1)) - zbands[:,np.newaxis,np.newaxis]), 0) < bandwid/2.
    ct2 = (X >= x_right) & (X < x_sink)
    ct = ct1 & ct2
    D[ct] = 0.8*n_crit
    
    # Left sink
    ct = (X > x_sink)
    D[ct] = 0.9*n_crit
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\slabcomb.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)

def slabpinch():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    
    x_crit = -17
    x_right = -14
    scale = 2.8
    n_crit = critDens(800)
    ct = (X < x_right) & (np.abs(Z) < 5.)
    D[ct] = n_crit * np.exp(-(-(X[ct] - x_crit))/scale)
    D[(D < 0.2*n_crit) & ct] = 0.2*n_crit
    

    # Middle bands
    zbands = np.array([-3,3]) # Location of bands
    bandwid = 0.2 # Width of bands
    
    x_sink = 3.0 # Location of the current sink
    ct1 = np.min(np.abs(np.tile(Z, (len(zbands),1,1)) - zbands[:,np.newaxis,np.newaxis]), 0) < bandwid/2.
    ct2 = (X >= x_right) & (X < x_sink)
    ct = ct1 & ct2
    D[ct] = 50*n_crit
    
    # Right sink
    ct = (X > x_sink)
    D[ct] = 50*n_crit
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\slabpinch.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)
    plt.axis("equal")
    
def stagesthree():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    
    x_wid = 4.0
    x_cent = -28
    ct1 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.)
    x_cent = -13.0
    ct2 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.)
    x_cent = +2.0
    ct3 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.)
    
    ct = (ct1 | ct2 | ct3) & (np.abs(Z) < 5.)
    D[ct] = 2.65e23 # Solid density target worth of electrons
    
    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stagesthree2.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)
    plt.axis("equal")

def stagesthreefuzz():
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    n_crit = critDens(800)
    
    x_wid = 1.0
    x_widb = 10.0
    x_cent = -28
    ct1 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.) & (np.abs(Z) < 5.)
    ct1b = (X < x_cent + x_widb/2.) & (X > x_cent - x_widb/2.) & (Z > 5.) & (Z < 7.)
    x_cent = -13.0
    ct2 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.) & (np.abs(Z) < 5.)
    ct2b = (X < x_cent + x_widb/2.) & (X > x_cent - x_widb/2.) & (Z < -5.) & (Z > -7.)
    x_cent = +2.0
    ct3 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.) & (np.abs(Z) < 5.)
    ct3b = (X < x_cent + x_widb/2.) & (X > x_cent - x_widb/2.) & (Z > 5.) & (Z < 7.)
    
    ct = (ct1 | ct2 | ct3) | (ct1b | ct2b | ct3b)
    D[ct] = 2.65e23 # Solid density target worth of electrons
    
    myscl = np.array([float(xmax - xmin)/nx, float(zmax - zmin)/nz])
    #D = gaussian_filter(D, radblur)
    Dblr1 = gaussian_filter(D, 0.1 / myscl)
    Dblr2 = gaussian_filter(D, 0.4 / myscl)
    D[np.abs(Z) < 6.0] = Dblr1[np.abs(Z) < 6.0]
    D[np.abs(Z) > 6.0] = Dblr2[np.abs(Z) > 6.0]


    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stagesthreefuzz.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)
    plt.axis("equal")

if __name__ == "__main__":
    zmin=-10
    zmax=10
    xmin=-30
    xmax=10
    npts=900
    nx = nz = npts
    xgv = np.linspace(xmin, xmax, nx)
    zgv = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xgv, zgv)
    
    D = np.zeros_like(X)
    n_crit = critDens(800)
    
    x_wid = 5.0
    x_widb = 12.0
    x_widc = 15.0
    x_cent = -28
    ct1 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.) & (np.abs(Z) < 3.5)
    ct1b = (X < x_cent + x_widb/2.) & (X > x_cent - x_widb/2.) & (Z > 3.5) & (Z < 4.5)
    #ct1c = (X < x_cent + x_widc/2.) & (X > x_cent - x_widc/2.) & (Z > 4.5) & (Z < 5)
    x_cent = -13.0
    ct2 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.) & (np.abs(Z) < 3.5)
    ct2b = (X < x_cent + x_widb/2.) & (X > x_cent - x_widb/2.) & (Z < -3.5) & (Z > -4.5)
    #ct2c = (X < x_cent + x_widc/2.) & (X > x_cent - x_widc/2.) & (Z < -4.5) & (Z > -5)
    x_cent = +2.0
    ct3 = (X < x_cent + x_wid/2.) & (X > x_cent - x_wid/2.) & (np.abs(Z) < 3.5)
    ct3b = (X < x_cent + x_widb/2.) & (X > x_cent - x_widb/2.) & (Z > 3.5) & (Z < 4.5)
    #ct3c = (X < x_cent + x_widc/2.) & (X > x_cent - x_widc/2.) & (Z > 4.5) & (Z < 5)
    
    ct = (ct1 | ct2 | ct3) | (ct1b | ct2b | ct3b)# | (ct1c | ct2c | ct3c) 
    D[ct] = 2.65e23 # Solid density target worth of electrons
    
    myscl = np.array([float(xmax - xmin)/nx, float(zmax - zmin)/nz])
    #D = gaussian_filter(D, radblur)
    Dblr1 = gaussian_filter(D, 0.1 / myscl)
    Dblr2 = gaussian_filter(D, 0.2 / myscl)
    D[np.abs(Z) < 4.0] = Dblr1[np.abs(Z) < 4.0]
    D[np.abs(Z) > 4.0] = Dblr2[np.abs(Z) > 4.0]


    fname = r"C:\Users\Scott\Documents\LSP\LSPwrite submissions\curtest1\stagesthreefuzz.dat"
    write2DLSP(fname, xgv, zgv, D/2.65)
    clabel = "Density"
    fig = plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    im = ax.pcolorfast((xmin, xmax), (zmin, zmax), D, vmin=0, vmax=n_crit, cmap='viridis')
    cbar = fig.colorbar(im, label=clabel)
    plt.axis("equal")
    

