#!/usr/bin/env python2
'''
Severe hack to get a line through the center of a sclr.p4 file.

Usage:
  sclrcenter.py [options] (--X|--Y|--Z) <input> (<var> <output>)...

Options:
  --help -h                Show this help.
  --verbose -v             Turn on verbosity.
  --X    -x                Use X line.
  --Y    -y                Use Y line.
  --Z    -z                Use Z line.
  --small=S -s S           Define the tolerance for what is close to the center. [default: 1.0e-6]
  --half                   If you KNOW you will get double, half it.
'''
import lspreader as rd;
import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.interpolate as interp;
from docopt import docopt;

def mkprint(name):
    def _print(ss):
        print("{}: {}".format(name,s));
    return _print;

def main():
    opts=docopt(__doc__,help=True);
    name = opts['<input>'];
    var = opts['<var>'];
    verbose = opts['--verbose'];
    sml = float(opts['--small']);
    logprint = mkprint(name) if verbose else lambda s: None;
    if opts['--X']:
        a,b = 'y','z';
    elif opts['--Y']:
        a,b = 'x','z';
    else:
        a,b = 'x','y';
    with rd.LspOutput(name, prefix=name) as f:#verbose=verbose
        d=f.get_data(var=var);
    amid = (d[a].max()-d[a].min())/2.0 + d[a].min();
    bmid = (d[b].max()-d[b].min())/2.0 + d[b].min();
    logprint("begin world's largest filter");
    good = (np.abs(d[a] - amid) < sml) & (np.abs(d[b] - bmid) < sml);
    logprint("done.");
    logprint("got {} matches.".format(good.astype(int).sum()));
    for var,outname in zip(opts['<var>'],opts['<output>']):
        logprint("outputting {}...".format(var));
        out = d[var][good];
        if opts['--half']: out = out[:out.size/2];
        out= {
            's':np.array([[out]]).T,
            'x':(d['x'].min(),d['x'].max()),
            'y':(d['y'].min(),d['y'].max()),
            'z':(d['z'].min(),d['z'].max())
        };
        with open(outname, "wb") as f:
            pickle.dump(out,f,2);
    pass;
if __name__=='__main__':
    main();
