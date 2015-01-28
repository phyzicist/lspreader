#!/usr/bin/env python2

import numpy as np;

def gendat(f,lims,filename=None):
    '''Generate a text file for custom functions for
       use with LSP. This is meant to be read as from
       a module.
       
       Parameters
       ==========
       f : A function. For now, this works for 2D data,
          but stay tuned, we'll make this work for
          arbitrary functions.
       lims : A list containing the limits of the sim-
          ulations which should be triples of min and
          max and steps.
       filename : Name of the file to output. If this
          is None, return a string. If this is a file
          object, append to the open file.'''
    if len(lims) != 2:
        s="I need to implement other dimensions. Be patient!"
        raise NotImplementedError(s);
    f = lambda x: np.linspace(x[0],x[1],x[2]);
    cords = map(f,lims);
    
    

    
