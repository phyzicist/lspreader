#!/usr/bin/env python2
'''
Convert an lsp fields or scalars into a 3d numpy array.

Usage:
  p4toS.py [options] <var> <randompickle> <input> <output>

Options:
  -h --help         Show this help.
  -v --verbose      Turn on verbosity.
  --pool_size=PS    Choose a pool size PS.
'''

import lspreader as rd;
import cPickle;
import numpy as np;
import scipy.interpolate as interpol;
from docopt import docopt;
import random;

def logprint(s):
    global verbose;
    if verbose:
        print("{}: {}".format(name,s));
    pass;

def interpolate_scalar(x,y,z,s,res=100j):
    '''Interpolates the scalar s on the grid x,y,z
       with the resolution res.'''
    Z,Y,X = np.mgrid[ min(z) : max(z) : res,
                      min(y) : max(y) : res,
                      min(x) : max(x) : res];
    S = interpol.griddata((x,y,z),s,(X,Y,Z));
    return S;

def main():
    global name,verbose;
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
    logprint('reading in {}'.format(inname));
    with rd.LspOutput(inname) as f:
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
    S = interpolate_scalar(x,y,z,s);
    del x,y,z,s,d;
    logprint("dumping");
    with open(outname,"wb") as f:
        cPickle.dump(S,f,2);
    pass;

if __name__ == '__main__':
    main();
