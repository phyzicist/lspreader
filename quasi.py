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
  --L=L -l L                  Give lambda for writing the critical density line. [default: 800e-9]
'''
import numpy as np;
import matplotlib;
import cPickle as pickle;
from matplotlib import colors;
from docopt import docopt;
from misc import read

e0=8.854187817e-12
m_e=9.10938291e-31
q=1.60217657e-19
c=2.99792458e8

def main():
    opts = docopt(__doc__,help=True);
    if opts['--output']: matplotlib.use('Agg');
    import matplotlib.pyplot as plt;
    if opts['--X']:
        get = lambda d: d[: ,d.shape[1]/2, d.shape[2]/2];
        label = 'x';
    elif opts['--Y']:
        get = lambda d: d[d.shape[0]/2, :, d.shape[2]/2];
        label = 'y';
    elif opts['--Z']:
        get = lambda d: d[d.shape[0]/2, d.shape[1]/2, :];
        label = 'z';
    Z = map(float, opts['<Z>']);
    names = zip(opts['<input>'],Z);
    #getting electrons first, using it to derive dimensions.
    e = read(opts['<electron-input>'],dumpfull=True);
    dim = e[label];
    e = get(e['s']);
    x = np.linspace(dim[0], dim[1], e.size);
    #setting up output
    nz = np.zeros(e.size);
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
    ncrit = e0*m_e/q**2*(2*np.pi*c/float(opts['--L']))**2;
    ncrit = np.log10(ncrit/1e6+0.1);
    plt.annotate('Critical Density', xy=((x[-1]-x[0])*0.1, ncrit+0.2 ));
    plt.axhline(ncrit,linestyle=':',c="gray");
    plt.plot(x, np.log10(e+0.1), label="$n_{ele}$",linewidth=3,linestyle=":",c="black");
    plt.plot(x, np.log10(nz+0.1), label=r"$n_{ion}\times Z$",linewidth=1);
    plt.ylabel("log$_{10}$ of number density (cm$^{-3}$)");
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
