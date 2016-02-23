# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:34:58 2016

This should cover the step between flash
Convert flash 3D simulation density file into a 2D LSP .dat
Save some plots

@author: Scott
"""


import yt
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate

# NumPy arrays
def findNearest(array, value):
    """
    Get the index for the nearest value to value
    Copied from: http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    Inputs:
        array: 1D array to search for value
        value: number to search for in array
    Outputs:
        nearest value in array, index of that value
    """
    idx = np.abs(array-value).argmin()
    return array[idx], idx
    
# YT slicing
def yslice(ds, fld = 'dens', lev = 5, y_um = 0.0):
    """ Extract a slice (fixed Y value) from 3D FLASH yt dataset.
    
    Note: I had to make this because Yt's function was putting vertical streaks in the data.
    As an added bonus, this one does some smoothing at the AMR level. Nice!
    
    Written by Scott Feister 2016-02-22.
    
    Inputs:
        ds: YT dataset. Loaded via ds = yt.load("myflash_hdf5_0000")
        fld: The field to look at, e.g. density
        lev: The refinement level up to which to extract. Bigger than max is fine, though array gets big.
        y_um: The y value (in microns) at which to extract the slice.
    Outputs:
        myslc: 2D (X by Z) numpy array containing extracted data
        xgv: 1D array, X values for the data
        zgv: 1D array, Z values for the data
    """
    
    px = ds.domain_width/ds.domain_dimensions # 1 (level 0) pixel, in unit code_length
    le = ds.domain_left_edge + 3 * px # Leave a 3 (level 0) pixel buffer to allow call to ds.smoothed_grid without yt error    
    re = ds.domain_right_edge - 3 * px # Leave a 3 (level 0) pixel buffer to allow call to ds.smoothed_grid without yt error    
    dims = np.array((re - le)/px) * ds.refine_by**lev # Specifies that an N x M x O pixel buffer will be made

    # Adjust such that we only take 2 values along Y dimension (N x 2 x O buffer in place of N x M x O; much smaller!)
    le[1] = ds.arr([y_um], 'um')[0].in_units('code_length') # We want the left edge to start at specified y value. Convert the input (in microns) to code_length units and inject into array. There is probably a better function, but hey, it works.
    dims[1] = 2 # Don't take all values along y dimension (memory overflow); take only 2. If we set this as 1, it for some reason will give back a blank.
    
    # Extract the smoothed grid data, if possible. Otherwise, extract raw grid data.
    try:
        gridded_data = ds.smoothed_covering_grid(level=lev, left_edge=le, dims=dims, fields=[fld])
    except:
        print "Error with smoothing function. Using non-smoothed version."
        gridded_data = ds.covering_grid(level=lev, left_edge=le, dims=dims, fields=[fld])
    
    # Do NumPy tricks to squeeze into a 2D array of the format I like ( X dims, Z dims ) long
    myslc = np.array(gridded_data[fld][:, 0, :]).swapaxes(0,1)
    
    # Specify the 1D grid vectors in X and Z for plotting, e.g. matplotlib.pyplot.pcolor(xgv, zgv, myslc)
    xgv = np.linspace(le[0], re[0], myslc.shape[1])*1e4
    zgv = np.linspace(le[2], re[2], myslc.shape[0])*1e4
    
    return myslc, xgv, zgv
    
def abelFLASH3D(ds):
    pass

def mypcolor(C, xgv, zgv, fig = None, cmin = 0,  cmax = None, title='', tstring = '', clabel = '', fld_id = '', sticker ='', rfooter = '', cmap='viridis'):
    """ A custom-tailored wrapper for pcolorfast(). Somewhat general, meant for any 2D FLASH colorplot.
    Inputs:
        C: 2D NumPy array, the data to be visualized
        xgv: x dimension grid vector (1D), in microns
        zgv: z dimension grid vector (1D), in microns
        ...
    """    
    # xgv, ygv should be in microns
    
    if not cmax:
        cmax = np.max(C)
    if not fig:
        fig = plt.figure()
    plt.clf() # Clear the figure
    ax = plt.subplot(111)
    xr = [xgv[0],xgv[-1]] # min and max of xgv, needed for pcolorfast
    zr = [zgv[0],zgv[-1]]
    plt.axis('equal')
    im = ax.pcolorfast(xr, zr, C, cmap=cmap)
    ax.set_xlabel(r'X ($\mu m$)')
    ax.set_ylabel(r'Z ($\mu m$)')
    ax.set_xlim(xr)
    ax.set_ylim(zr)
    ax.set_title(title, fontsize=20)
    ax.text(np.min(xgv) + 1, np.max(zgv) - 3, tstring, fontsize=24, color='white')
    ax.text(np.max(xgv) - 6, np.min(zgv) + 2, fld_id, fontsize=44, color='white')
    ax.text(np.min(xgv) + 1, np.min(zgv) + 2, sticker, fontsize=44, color='white')
    cbar = fig.colorbar(im, label=clabel)
    im.set_clim(vmin=cmin, vmax=cmax)
    fig.text(0.99, 0.01, rfooter, horizontalalignment='right')
    #sp.addCrit(ax, C, zgv, rgv)
    return fig

def writeLSP(fname, dens, dens_zgv, dens_xgv, pc_xdims, pc_zdims, nxpoints = 400, nzpoints = 403, xshift=-999):
    """ Write an LSP 2D .dat output from XYZ 3D FLASH inputs 
    """
    # Convert dimensions back into LSP units    
    (xmin, xmax) = np.array(pc_xdims)*1e-4 # LSP X dimensions in cm (laser prop direction)
    (zmin, zmax) = np.array(pc_zdims)*1e-4 # LSP Z dimensions in cm
    dens_zgv = dens_zgv*1e-4 # LSP Z dimensions in cm
    dens_xgv = dens_xgv*1e-4 # LSP X dimensions in cm (laser prop direction)
        
    xgv = np.linspace(xmin, xmax, nxpoints)
    zgv = np.linspace(zmin, zmax, nzpoints)
    X, Z = np.meshgrid(xgv, zgv)

    myf = interpolate.RectBivariateSpline(dens_zgv, dens_xgv, dens)
    D = myf(np.abs(Z), X, grid=False)

    with open(fname, 'w') as f:
        f.write('# FLASH-generated LSP dat file (function type 40). Corona 1.3 microns; defocused. Original target surface at X = ' + str(xshift) + ' microns. \n')
        f.write('# Dimensions X, Z, Density. See page 147 of 200 in LSP manual, or search "type 40", for specification format.\n')
        f.write(str(int(nxpoints)) + ' ' + str(int(nzpoints)) + '\n')
        np.savetxt(f, xgv[None], delimiter = ' ') # xgv[None] changes the dimensions of xgv from (100,) to (1, 100); this lets us save it as a row rather than column
        np.savetxt(f, zgv[None], delimiter = ' ')
        np.savetxt(f, D, delimiter = ' ')

    return D, Z*1e4, X*1e4

def FLASHtoLSP(fn_h5, fn_dat, focx_flash, pc_xdims, pc_zdims, plotdir, shortname = '',  ncrit = 1.74e21, xcrit_lsp = -6.1):
    """ Custom convert 3D FLASH HDF5 data file into an LSP .dat, specifically for Scott's parameter scan
    Inputs:    
        fn_h5: string, full filename of input FLASH 3D HDF5
        fn_dat: string, full filename for output LSP 2D .dat
        focx_flash: number, laser focal point X value in FLASH simulation
        ncrit: number, Electron critical density in elec/cc
        xcrit_lsp: number, Desired X value (microns) for 6-ionized critical density, for LSP
        pc_xdims: 1D, two-element array. X limits of LSP particle creation space
        pc_zdims: 1D, two-element array. Z limits of LSP particle creation space
        plotdir: string, full folder path for output plots of FLASH to LSP conversion
    Outputs:
        focx_lsp: number, laser focal point X value to set in LSP simulation
        Writes a LSP 2D .dat to fn_dat
    """
    # Load in the FLASH data and extract density slice at Y=0
    ds = yt.load(fn_h5)
    dens, xgv, zgv = yslice(ds, fld = 'dens') # Density is in gram/cc. xgv and zgv refer to FLASH dimensions.
    
    # Extract the simulation timestamp
    time_ns = float(ds.current_time.in_units('ns')) # Current simulation time
    print "Timestamp of FLASH HDF5 file:", time_ns, "ns"
    
    # Convert density to number/cc (singly and fully-ionized versions)    
    dens = dens * 1.0e23 # Initial electron density (1 electron per Oxygen)
    fulldens = 6 * dens # Fully ionized electron density (6 electrons per Oxygen)

    # Extract a central lineout of the fully-ionized electron density (near Z=0)
    cidx = findNearest(zgv, 0)[1] # Index of (or nearest to) Z = 0
    fdline = fulldens[cidx, :] # X=0 lineout of the fully-ionized density profile

    
    ixcrit = np.argmax(fdline > ncrit) # Array index of the critical density, when full-ionized
    xcrit_flash = xgv[ixcrit] # X value of six-ionized critical density, in  FLASH coordinates
    print xcrit_flash
    
    xshift = xcrit_lsp - xcrit_flash
    print "X shift", xshift
    xgv_new = xgv + xshift

    # Calculate how the focal depth should now shift
    #print "Focal depth (um) past target original surface:", focx 
    focx_lsp = focx_flash + xshift
    print "Focal depth X (um) to code into this LSP sim", focx_lsp
    
    # Convert to an LSP .dat
    D, Z, X = writeLSP(fn_dat, dens, zgv, xgv_new, pc_xdims, pc_zdims, xshift=xshift)
    
    ## Make some plots, save to PNG
    # Subcritical FLASH 2D plot
    fig = plt.figure(1)
    mypcolor(fulldens, xgv, zgv, fig, cmax = ncrit, title='FLASH density * six-ionized', clabel = 'Density (elec/cc)', rfooter = shortname)
    fig.savefig(os.path.join(plotdir, "FLASH subcritical.png"))    
    
    # Overdense FLASH 2D plot
    fig = plt.figure(2)
    mypcolor(dens, xgv, zgv, fig, cmax = 2.0e23, title='FLASH density * single-ionized', clabel = 'Density (elec/cc)', rfooter = shortname)
    fig.savefig(os.path.join(plotdir, "FLASH overdense.png"))    

    # LSP central lineout
    fig = plt.figure(4)
    plt.clf()
    plt.plot(xgv_new, fdline)
    plt.ylim(0, ncrit)
    plt.title("FLASH density central lineout")
    plt.xlabel("Laser direction (LSP X coordinate) (microns)")
    plt.ylabel("Density, six-ionized (elec/cc)")
    plt.xlim(xgv_new[0], xcrit_lsp)
    fig.text(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, "FLASH subcritical lineout.png"))    

    # LSP subcritical 2D plot
    fig = plt.figure(5)
    mypcolor(D, X[0,:], Z[:,0], fig, cmax = 1.7e21, title="LSP input", clabel = 'Density (elec/cc)', rfooter = shortname)
    fig.savefig(os.path.join(plotdir, "FLASH-generated LSP input, subcritical.png"))    

    plt.close('all') # Close the figures
    
    return focx_lsp
#    
#if __name__ == "__main__":
#    ## EXAMPLE CODE ONLY. WILL NOT ACTUALLY SUCCESSFULLY RUN.
#    shortname = r"3Dtest10_f15c13"
#    runroot = 'C:\Users\Scott\Documents\temp\FlashRams'
#    plotroot = 'C:\Users\Scott\Documents\temp\FlashRams\Myplots'
#
#    rundir = sf.subdir(runroot, shortname)
#    plotdir = sf.subdir(plotroot, shortname)
#    
#    # Get the Checkpoint file name
#    fn_h5 = sf.getfns(rundir, prefix = "lasslab_hdf5_chk_")[-1]
#    
#    # Check that we're loading a checkpoint from after the FLASH sim has actually run
#    if fn_h5[-4:] == "0000":
#        raise Exception("FLASH checkpoint 0000 is largest checkpoint; FLASH sim hasn't completed? Aborting FLASH to LSP conversion.")
#    
#    fn_dat = os.path.join(rundir, 'watercolumn.dat')
#
#    focx = 15.0 # Focal depth, in coordinates where X = 0 at the unperturbed water column surface.
#
#    focx_flash = focx - 15 # Focal depth, X in FLASH coordinates. Account for the fact that FLASH X = 0 is center, not surface, of target. Change such that surface of target (15 micron column) is now centered.
#
#    pc_xdims = [-30, 10]
#    pc_zdims = [-15, 15]
#
#    focx_lsp = FLASHtoLSP(fn_h5, fn_dat, focx_flash, pc_xdims, pc_zdims, plotdir, ncrit = 1.74e21, xcrit_lsp = -6.1)
