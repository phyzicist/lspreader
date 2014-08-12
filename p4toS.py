#!/usr/bin/env python2
'''
Convert an lsp fields or scalars into a 3d numpy array.

Usage:
  p4toS.py [options] <var> <randompickle> <input> <output>

Options:
  -h --help         Show this help.
  -v --verbose      Turn on verbosity.
  --pool_size=PS    Choose a pool size PS.
  -x --X            Use X in interpolation
  -y --Y            Use Y in interpolation
  -z --Z            Use Z in interpolation
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

def interpolate_scalar_3d(x,y,z,s,res=100j):
    '''Interpolates the scalar s on the grid x,y,z
       with the resolution res.'''
    Z,Y,X = np.mgrid[ min(z) : max(z) : res,
                      min(y) : max(y) : res,
                      min(x) : max(x) : res];
    S = interpol.griddata((x,y,z),s,(X,Y,Z));
    return S;

def interpolate_scalar_2d(x,y,s,res=100j):
    '''Interpolates the scalar s on the grid x,y
       with the resolution res.'''
    Y,X = np.mgrid[min(y) : max(y) : res,
                   min(x) : max(x) : res];
    S = interpol.griddata((x,y),s,(X,Y));
    return S;

def interpolate_scalar_1d(x,s,res=100j):
    '''Interpolates the scalar s on the grid x
       with the resolution res.'''
    X = np.mgrid[min(x) : max(x) : res];
    S = interpol.griddata((x,),s,(X,));
    return S;

def process_args():
    opts=docopt(__doc__,help=True);
    var=opts['<var>'];
    rndname=opts['<randompickle>'];
    name = inname = opts['<input>'];
    outname = opts['<output>'];
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
    return var,rndname,name,verbose,outname,pool_size,useX,useY,useZ;

def main():
    global name,verbose;
    var,rndname,name,verbose,outname,pool_size,useX,useY,useZ=process_args();
    logprint('reading in {}'.format(name));
    with rd.LspOutput(name, verbose=verbose, prefix=name) as f:
        d = f.get_data(var=[var],pool_size=pool_size,lazy=False);
    logprint("preparing to prune, reading in {}".format(rndname));
    with open(rndname,'rb') as f:
        rnd = cPickle.load(f);
    logprint("pruning...");
    for k in d:
        d[k] = [d[k][i] for i in rnd];
    logprint('making arrays for interpolating scalar field {}'.format(var));
    x=np.array(d['x']);
    y=np.array(d['y']);
    z=np.array(d['z']);
    s=np.array(d[var]);
    logprint('interpolating');
    if useX and useY and useZ:
        S = interpolate_scalar_3d(x,y,z,s);
    elif useX and useY:
        S = interpolate_scalar_2d(x,y,s);
    elif useX and useZ:
        S = interpolate_scalar_2d(x,z,s);
    elif useY and useZ:
        S = interpolate_scalar_2d(y,z,s);
    elif useX:
        S = interpolate_scalar_1d(x,s);
    elif useY:
        S = interpolate_scalar_1d(y,s);
    elif useZ:
        S = interpolate_scalar_1d(z,s);
    del x,y,z,s,d;
    logprint("dumping");
    with open(outname,"wb") as f:
        cPickle.dump(S,f,2);
    pass;

if __name__ == '__main__':
    main();
