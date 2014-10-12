#!/usr/bin/env python2
'''
Convert an lsp fields or scalars into a 3d numpy array.

Usage:
  sclr.py [options] <input> (<var> <output>)...

Options:
  -h --help               Show this help.
  -v --verbose            Turn on verbosity.
  -x --X                  Use X in interpolation.
  -y --Y                  Use Y in interpolation.
  -z --Z                  Use Z in interpolation.
  --xres=XRES             Set the resolution along the x direction [default: 100].
  --yres=XRES             Set the resolution along the y direction [default: 100].
  --zres=XRES             Set the resolution along the z direction [default: 100].
'''

import lspreader as rd;
import cPickle;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.interpolate as interpol;
from docopt import docopt;
import random;
from misc import conv;

def logprint(s):
    global verbose;
    if verbose:
        print("{}: {}".format(name,s));
    pass;

def histogram_scalar_3d(x,y,z,s,
                        xres=100,yres=100,zres=100):
    '''Histograms the scalar s on the grid x,y,z
       with the resolutions *res along each axis.'''
    xbins = np.linspace(min(x),max(x),xres+1);
    ybins = np.linspace(min(y),max(y),yres+1);
    zbins = np.linspace(min(z),max(z),zres+1);
    H,_= np.histogramdd([x,y,z],bins=[xbins,ybins,zbins],weights=s);
    #divide to get the average
    norm=float(len(s))/float(xres*yres*zres);
    logprint('we have a norm of {}'.format(norm));
    return H/norm;

def histogram_scalar_2d(x,y,s,
                        xres=100,yres=100):
    '''Histograms the scalar s on the grid x,y
       with the resolutions *res along each axis.'''
    xbins = np.linspace(min(x),max(x),xres+1);
    ybins = np.linspace(min(y),max(y),yres+1);
    
    H,_,_ = np.histogram2d(x,y,bins=(xbins,ybins),weights=s);
    norm=float(len(s))/float(xres*yres);
    logprint('we have a norm of {}'.format(norm));
    return H/norm;

def histogram_scalar_1d(x,s,res=100):
    '''Histograms the scalar s on the grid x
       with the resolution res.'''
    xbins = np.linspace(min(x),max(x),res+1);
    H,_ = np.histogram(x,bins=xbins,weights=s);
    norm=float(len(s))/float(res);
    logprint('we have a norm of {}'.format(norm));
    return H/norm;

def main():
    global name,verbose;
    #procesing arguments
    opts=docopt(__doc__,help=True);
    var=opts['<var>'];
    outnames= opts['<output>'];
    vopairs = zip(var,outnames);
    name = opts['<input>'];

    verbose = opts['--verbose'];
    use = [];
    res = [];
    if opts['--X']:
        use.append('x');
        res.append(int(opts['--xres']));
    if opts['--Y']:
        use.append('y');
        res.append(int(opts['--yres']));
    if opts['--Z']:
        use.append('z');
        res.append(int(opts['--zres']));
    if use == []:
        use = ['x','y','z'];
        res = map(lambda k: int(opts[k]),['--xres','--yres','--zres']);
    # A couple of things to note; written in this way, whatever
    # this list (and thus, what is read) becomes, it is ordered
    # alphabetically. This is important, as this determines what
    # each resulting row and column and breadth in the output
    # array corresponds to from the actual simulation.
    #
    # It is probably worth mentioning that the xz in simulation
    # axes will be [0,1] in numpy axes, that is, it will be left-handed.
    # Using xz leads to this anyway, but it's worth reminding the reader.
    logprint('reading in {}'.format(name));
    with rd.LspOutput(name, verbose=verbose, prefix=name) as f:
        d = f.get_data(var=var);
    for v,outname in vopairs:
        logprint('histogramming');
        if len(use) == 3:
            S = histogram_scalar_3d(d['x'],d['y'],d['z'],s,
                                    xres=res[0],yres=res[1],zres=res[2]);
        elif len(use) == 2:
            S = histogram_scalar_2d(d[use[0]],d[use[1]],s,
                                    xres=res[0],yres=res[1]);
        else:
            S = histogram_scalar_1d(d[use[0]],s,res=res[0]);
        logprint("dumping");
        with open(outname,"wb") as f:
            cPickle.dump(S,f,2);
    pass;

if __name__ == '__main__':
    main();
