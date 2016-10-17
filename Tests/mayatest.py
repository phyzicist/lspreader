#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
mayatest.py: Test of the MayaVi package

Created by Scott Feister on Sun Oct 16 20:26:47 2016
"""
import mayavi
import mayavi.mlab as mlab
import numpy as np
mlab.options.offscreen = True
from scipy.interpolate import RectBivariateSpline

def example1():
    """ A simple 3D plot, saved to a png """
    fnout = r"C:\Users\Scott\Documents\temp\oct2016\OUTS\mlab tests\test.png"
    mlab.clf() # Clear current figure
    
    X, Y = np.mgrid[-7.:7.05:0.1, -5.:5.05:0.05]
    C = np.sin(X + Y) + np.sin(2 * X - Y) + np.cos(3 * X + 4 * Y)    
    s = mlab.surf(X,Y,C,colormap='gist_earth', warp_scale=0.5)
    
    xpts = np.array([1])
    ypts = np.array([3])
    zpts = np.array([2])
    mlab.axes() # Turns on the X/Y/Z axes
    mlab.points3d(xpts,ypts,zpts)
    
    mlab.savefig(fnout)
    mlab.close()
    return

def example2():
    """ More advanced example of plotting with various spot sizes """
    fnout = r"C:\Users\Scott\Documents\temp\oct2016\OUTS\mlab tests\test2.png"
    mlab.clf() # Clear current figure
        
    # Generate fake scalar field data
    X, Y = np.mgrid[-1.5:1.55:0.02, -1.5:1.55:0.02]
    C = np.sin(X + Y) + np.sin(2 * X - Y) + np.cos(3 * X + 4 * Y)    
    
    # Generate fake scatter points data
    npts = 100 # Number of points
    xpts = 2*(0.5 - np.random.rand(100))
    ypts = 2*(0.5 - np.random.rand(100))
    spts = np.random.rand(100) # Sizes; Not sure if this is area or what
    cpts = np.sqrt(xpts**2 + ypts**2) # Colormap values


    # Prep plot
    ws = 0.5 # Warp_scale
    s = mlab.surf(X,Y,C,colormap='gist_earth', warp_scale=ws)
    mlab.axes() # Turns on the X/Y/Z axes

    # Interpolate to get dots down onto point
    rbv = RectBivariateSpline(X[:,0], Y[0,:], C)
    zpts = rbv.ev(xpts, ypts)*ws # The height will always match the field

    # Add scatter points
    pt3d = mlab.points3d(xpts,ypts,zpts, cpts, scale_factor=0.1, colormap='Purples') # Assign colors
    pt3d.glyph.scale_mode = 'scale_by_vector' # Assign sizes, step 1
    pt3d.mlab_source.dataset.point_data.vectors = np.tile(spts, (3,1)).T # Assign sizes, step 2

    # Save figure and close
    mlab.savefig(fnout)
    mlab.close()
    return
    
if __name__ == "__main__":
    example2()
