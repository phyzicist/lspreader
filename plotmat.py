#!/usr/bin/env python2
import numpy as np;
import matplotlib.pyplot as plt;
import cPickle;
import math;
from docopt import docopt;

usage='''Plot cutplane at center

Usage:
  ./plotmat.py <infile> <outfile> <time>
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
    S = read_file(opts['<infile>']);
    S = zero_nan(S);
    X,Y = np.mgrid[-20:20:100j,-30:5:100j];
    plt.pcolormesh(X,Y,S[:,50,:]);
    plt.xlabel('x ($\mu m$)');
    plt.ylabel('z ($\mu m$)');
    plt.annotate('t = {} fs'.format(opts['<time>']),(-19,3),color='w');
    c=plt.colorbar();
    c.set_label('log10 of density');
    plt.savefig(opts['<outfile>']);
pass;
if __name__ == '__main__':
    main();
