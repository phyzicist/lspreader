#!/usr/bin/env python2
'''
Convert an lsp fields or scalars into a 3d numpy array.

Usage:
  p4toH.py [options] <input> (<var> <output>)...

Options:
  -h --help               Show this help.
  -v --verbose            Turn on verbosity.
  -p PS --pool_size=PS    Choose a pool size PS.
  -x --X                  Use X in interpolation
  -y --Y                  Use Y in interpolation
  -z --Z                  Use Z in interpolation
'''

import lspreader as rd;
import cPickle;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.interpolate as interpol;
from docopt import docopt;
import random;

def logprint(s):
    global verbose;
    if verbose:
        print("{}: {}".format(name,s));
    pass;

def histogram_scalar_3d(x,y,z,s,res=100):
    '''Histograms the scalar s on the grid x,y,z
       with the resolution res.'''
    #Z,Y,X = np.mgrid[ min(z) : max(z) : res*1j,
    #                  min(y) : max(y) : res*1j,
    #                  min(x) : max(x) : res*1j];
    xbins = np.linspace(min(x),max(x),res+1);
    ybins = np.linspace(min(y),max(y),res+1);
    zbins = np.linspace(min(z),max(z),res+1);
    H,_= np.histogramdd([x,y,z],bins=[xbins,ybins,zbins],weights=s);
    return H;

def histogram_scalar_2d(x,y,s,res=100):
    '''Histograms the scalar s on the grid x,y
       with the resolution res.'''
    xbins = np.linspace(min(x),max(x),res+1);
    ybins = np.linspace(min(y),max(y),res+1);
    
    H,_,_ = np.histogram2d(x,y,bins=(xbins,ybins),weights=s);
    return H;

def histogram_scalar_1d(x,s,res=100):
    '''Histograms the scalar s on the grid x
       with the resolution res.'''
    xbins = np.linspace(min(x),max(x),res+1);
    H,_ = np.histogram(x,bins=xbins,weights=s);
    return H;

def main():
    global name,verbose;
    #procesing arguments
    opts=docopt(__doc__,help=True);
    var=opts['<var>'];
    outnames= opts['<output>'];
    vopairs = zip(var,outnames);
    name = inname = opts['<input>'];

    verbose = opts['--verbose'];
    if opts['--pool_size']:
        pool_size=int(opts['--pool_size']);
    else:
        pool_size=24;
    if opts['--X']:
        useX=True;
    else:
        useX=False;
    if opts['--Y']:
        useY=True;
    else:
        useY=False;
    if opts['--Z']:
        useZ=True;
    else:
        useZ=False;
    if not useX and  not useY and not useZ:
        useX=useY=useZ=True;
    logprint('reading in {}'.format(name));
    with rd.LspOutput(name, verbose=verbose, prefix=name) as f:
        d = f.get_data(var=var,pool_size=pool_size);
    for v,outname in vopairs:
        logprint('making arrays for interpolating scalar field {}'.format(v));
        x=np.array(d['x']);
        y=np.array(d['y']);
        z=np.array(d['z']);
        s=np.array(d[v]);
        logprint('histogramming');
        if useX and useY and useZ:
            S = histogram_scalar_3d(x,y,z,s);
        elif useX and useY:
            S = histogram_scalar_2d(x,y,s);
        elif useX and useZ:
            S = histogram_scalar_2d(x,z,s);
        elif useY and useZ:
            S = histogram_scalar_2d(y,z,s);
        elif useX:
            S = histogram_scalar_1d(x,s);
        elif useY:
            S = histogram_scalar_1d(y,s);
        elif useZ:
            S = interpolate_scalar_1d(z,s);
        logprint("dumping");
        with open(outname,"wb") as f:
            cPickle.dump(S,f,2);
    pass;

if __name__ == '__main__':
    main();
