#!/usr/bin/env python2
'''Plot a cutplane of a 3D scalar file.

Usage:
  ./plotmat.py [options] [(--X | --Y | --Z)] [(--half|--index=INDEX)] <infile> <outfile> <time>

Options:
  --min=MIN -n MIN            Plot with a minimum MIN.
  --max=MAX -x MAX            Plot with a minimum MAX.
  --X                         Plot with X set to the index. Default.
  --Y                         Plot with Y set to the index.
  --Z                         Plot with Z set to the index.
  --index=INDEX -i INDEX      The index that defines the cutplane.
  --half                      Just cut it in half.
  --T                         Transpose the SP.
'''


import numpy as np;
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt;
import cPickle;
import math;
from docopt import docopt;


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
    opts=docopt(__doc__,help=True);
    vmin = float(opts['--min']) if opts['--min'] else -1.0;
    vmax = float(opts['--max']) if opts['--max'] else 23.5;
    if not opts['--X'] and not opts['--Y'] and not opts['--Z']:
        opts['--X']=True;
    S = read_file(opts['<infile>']);
    S = zero_nan(S);
    S = np.log10(S+0.1);
    if opts['--index']:
        i = int(opts['--index']);
    elif opts['--half']:
        if opts['--X']:
            i = len(S[:,0,0])/2;
        elif opts['--Y']:
            i = len(S[0,:,0])/2;
        else:
            i = len(S[0,0,:])/2;
    else:
        i = 0;
    if opts['--X']:
        SP = S[i,:,:];
    elif opts['--Y']:
        SP = S[:,i,:];
    else:
        SP = S[:,:,i];
    if opts['--T']:
        SP = SP.T;
    X,Y = np.mgrid[-20:20:100j,-30:5:100j];
    plt.pcolormesh(X,Y,SP,vmin=vmin,vmax=vmax);
    plt.xlabel('z ($\mu m$)');
    plt.ylabel('x ($\mu m$)');
    plt.annotate('t = {} fs'.format(opts['<time>']),(-19,3),color='w');
    c=plt.colorbar();
    c.set_label('log10 of density');
    plt.savefig(opts['<outfile>']);
pass;
if __name__ == '__main__':
    main();
