#!/usr/bin/env python2
import numpy as np;
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt;
import cPickle;
import math;
from docopt import docopt;

usage='''Plot cutplane at center

Usage:
  ./plotmat.py [--min=MIN --max=MAX] <infile> <outfile> <time>
'''

def read_file(filename):
    print('loading file {}'.format(filename));
    with open(filename,'rb') as f:
        d=cPickle.load(f);
    return d;

def tozero(v):
    if math.isnan(v):
        return 0;
    else:
        return v;

def zero_nan(S):
    return np.array([[[tozero(k) for k in j] for j in i] for i in S]);

def main():
    opts=docopt(usage,help=True);
    if opts['--min']:
        vmin = float(opts['--min']);
    else:
        vmin = -1.0;
    if opts['--max']:
        vmax = float(opts['--max']);
    else:
        vmax = 23.5;
    S = read_file(opts['<infile>']);
    S = zero_nan(S);
    X,Y = np.mgrid[-20:20:100j,-30:5:100j];
    plt.pcolormesh(X,Y,S[:,50,:],vmin=vmin,vmax=vmax);
    plt.xlabel('z ($\mu m$)');
    plt.ylabel('x ($\mu m$)');
    plt.annotate('t = {} fs'.format(opts['<time>']),(-19,3),color='w');
    c=plt.colorbar();
    c.set_label('log10 of density');
    plt.savefig(opts['<outfile>']);
pass;
if __name__ == '__main__':
    main();
