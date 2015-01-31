#!/usr/bin/env python2
'''
Plot the n_ion*Z plot.

Usage:
  quasi.py [options] (--X|--Y|--Z) (<input> <Z>)...

Options:
  --title=TITLE -t Title      Set the title.
  --X -x                      Draw a line along X, half in the other dimensions.
  --Y -y                      Draw a line along Y, half in the other dimensions.
  --Z -z                      Draw a line along Z, half in the other dimensions.
  --restrict=R                Restrict to range R.
  --high-res -H               Output a high resolution plt.
  --ylim=YLIM                 Set ylim for videos.
  --factor=F -f F             Multiply histogram by F. [default: 1.0]
'''
import numpy as np;
import matplotlib.pyplot as plt;
import cPickle as pickle;
from matplotlib import colors;
from docopt import docopt;
from misc import read
def main():
    opts = docopt(__doc__,help=True);
    if opts['--X']:
        get = lambda d: d[: ,d.shape[1]/2, d.shape[2]/2];
        l = 'x';
    elif opts['--Y']:
        get = lambda d: d[d.shape[0]/2, :, d.shape[2]/2];
        l = 'y';
    elif opts['--Z']:
        get = lambda d: d[d.shape[0]/2, d.shape[1]/2, :];
        l = 'z';
    Z = map(float, opts['<Z>']);
    names = zip(opts['<input>'],Z);
    #handling first
    name,z = names[0]; names = names[1:];
    d = read(name,dumpfull=True);
    n  = get(d['s']);
    nz = n*z;
    x = np.linspace(d[l][0], d[l][1], n.size);
    #looping over files
    for name,z in zip(opts['<input>'],Z):
        curn  = get(read(name));
        n += curn
        nz+= curn*z;
    #restricting
    if opts['--restrict']:
        r = eval(opts['--restrict']);
        good = (x > r[0]) & (x < r[1]);
        x  =  x[good];
        n  =  n[good];
        nz = nz[good];
        print(len(x));
    #plotting
    plt.plot(x, np.log10(n+0.1), label="$n_{ion}$");
    plt.plot(x, np.log10(nz+0.1), label="$n_{ion}Z$");
    plt.xlim(x[0],x[-1]);
    if opts['--ylim']:
        ylim=eval(opts['--ylim']);
        plt.ylim(ylim);
    plt.legend(loc='upper left');
    plt.show();
    
if __name__ == '__main__':
    main();
