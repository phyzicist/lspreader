#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
georeader.py: Read geometry files from LSP outputs

Created by Scott Feister on Thu Oct 20 13:27:08 2016

Readers for grid.p4, ... written for LSP-10

CHANGELOG:
-
TODO:
* Check stuff
"""

from collections import OrderedDict
import numpy as np

# TODO: Create something a bit more organized as far as outputs go? Maybe it's ok.
def read_grid(fn):
    """ Read an (always-ASCII) 'grid.p4' LSP output file into a python dictionary
    fn: filename (e.g. fn = "/home/joeblow/myrun/grid.p4")
    
    """
    d = OrderedDict()
    with open(fn, 'r') as f:
        d['title'] = f.readline().strip()
        d['GEOMETRY'] = int(f.readline().strip())
        d['DIMENSION'] = int(f.readline().strip())
        d['x_units'] = f.readline().strip()
        d['y_units'] = f.readline().strip()
        d['z_units'] = f.readline().strip()
        d['ngrids'] = int(f.readline().strip()) # Number of grids
        
        grids = [OrderedDict()]*d['ngrids'] # Make a list of grids
        d['grids'] = grids
        for i in range(d['ngrids']): # For each grid
            grids[i]['index'] = int(f.readline().strip())
            
            nI = int(f.readline().strip())
            grids[i]['nI'] = nI
            grids[i]['xgv'] = np.zeros(nI)
            for j in range(nI):
                grids[i]['xgv'][j] = float(f.readline().strip())
                
            nJ= int(f.readline().strip())
            grids[i]['nJ'] = nJ
            grids[i]['ygv'] = np.zeros(nJ)
            for j in range(nJ):
                grids[i]['ygv'][j] = float(f.readline().strip())
                
            nK= int(f.readline().strip())
            grids[i]['nK'] = nK
            grids[i]['zgv'] = np.zeros(nK)
            for j in range(nK):
                grids[i]['zgv'][j] = float(f.readline().strip())
    return d

# TODO: Test on a more-than-one region
def read_regions(fn):
    """ Read an (always-ASCII) 'regions.p4' LSP output file into a python dictionary """
    d = OrderedDict()
    with open(fn, 'r') as f:
        d['title'] = f.readline().strip()
        d['GEOMETRY'] = int(f.readline().strip())
        d['DIMENSION'] = int(f.readline().strip())
        d['nregions'] = int(f.readline().strip()) # Number of regions
        regs = [OrderedDict()]*d['nregions'] # Make a list of regions
        d['regions'] = regs
        for i in range(d['nregions']): # For each region
            line = f.readline().strip()
            myarr = np.fromstring(line, dtype='int', sep=" ", count=3)
            (regs[i]['nI'], regs[i]['nJ'], regs[i]['nK']) = myarr
    
            line = f.readline().strip()
            myarr = np.fromstring(line, dtype='float', sep=" ", count=6)
            (regs[i]['xmin'], regs[i]['xmax'], regs[i]['ymin'], regs[i]['ymax'], regs[i]['zmin'], regs[i]['zmax']) = myarr
    return d

# TODO: Test on a non-zero volume
def read_volumes(fn):
    """ Read an (always-ASCII) 'volumes.p4' LSP output file into a python dictionary """
    d = OrderedDict()
    with open(fn, 'r') as f:
        d['title'] = f.readline().strip()
        d['nvolumes'] = int(f.readline().strip()) # Number of regions
        vols = [OrderedDict()]*d['nvolumes'] # Note: If no volumes, this returns an empty list, which is fine -- just FYI.
        d['volumes'] = vols
        for i in range(d['nvolumes']):
            line = f.readline().strip()
            myarr = np.fromstring(line, dtype='float', sep=" ", count=7)
            (_, vols[i]['x0'], vols[i]['x1'], vols[i]['y0'], vols[i]['y1'], vols[i]['z0'], vols[i]['z1']) = myarr
            myarr = np.fromstring(line, dtype='int', sep=" ", count=1)
            (vols[i]['type'],) = myarr
    return d
    
if __name__ == "__main__":
    fn_grid = r"C:\Users\Scott\Documents\temp\oct2016\test3D_0\grid.p4"
    d_grid = read_grid(fn_grid)
    
    fn_regs = r"C:\Users\Scott\Documents\temp\oct2016\test3D_0\regions.p4"
    d_regs = read_regions(fn_regs)

    fn_vols =  r"C:\Users\Scott\Documents\temp\oct2016\test3D_0\volumes.p4"
    d_vols = read_volumes(fn_vols)