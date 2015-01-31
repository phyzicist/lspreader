#!/usr/bin/env python2
'''
Plot the n_ion*Z plot.

Usage:
  quasi.py [options] (--X|--Y|--Z) <electron-input> (<input> <Z>)...

Options:
  --title=TITLE -t TITLE      Set the title.
  --output=OUT -o OUT         Output instead of display.
  --X -x                      Draw a line along X, half in the other dimensions.
  --Y -y                      Draw a line along Y, half in the other dimensions.
  --Z -z                      Draw a line along Z, half in the other dimensions.
  --restrict=R                Restrict to range R.
  --high-res -H               Output a high resolution plt.
  --ylim=YLIM                 Set ylim for videos.
  --factor=F -f F             Multiply histogram by F. [default: 1.0]
'''
import numpy as np;
import matplotlib;
import cPickle as pickle;
from matplotlib import colors;
from docopt import docopt;
from misc import read
def main():
    opts = docopt(__doc__,help=True);
    if opts['--output']: matplotlib.use('Agg');
    import matplotlib.pyplot as plt;
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
    #getting electrons first, using it to derive dimensions.
    e = read(opts['<electron-input>'],dumpfull=True);
    x = np.linspace(e[l][0], e[l][1], e.size);
    e = get(e['s']);
    #setting up output
    nz = zeros(e.size);
    #looping over other species
    for name,z in zip(opts['<input>'],Z):
        nz += get(read(name))*z;
    #restricting
    if opts['--restrict']:
        r = eval(opts['--restrict']);
        good = (x > r[0]) & (x < r[1]);
        x  =  x[good];
        e  =  e[good]
        nz = nz[good];
    #plotting
    plt.plot(x, e, label="$n_{ele}$");
    plt.plot(x, nz, label="$n_{ion}\timesZ$");
    plt.ylabel("log$_{10}$ of number density (");
    if opts['--title']: plt.title(opts['--title']);
    plt.xlim(x[0],x[-1]);
    if opts['--ylim']:
        ylim=eval(opts['--ylim']);
        plt.ylim(ylim);
    plt.legend(loc='upper left');
    if opts['--output']:
        plt.savefig(opts['--output']);
    else:
        plt.show();
    
if __name__ == '__main__':
    main();
