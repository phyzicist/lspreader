#!/usr/bin/env python2
'''Plot a cutplane of a 3D scalar file.

Usage:
  ./plotcut.py [options] [(--X | --Y | --Z)] [(--half|--index=INDEX)] <infile> [<outfile>]

Options:
  --min=MIN -n MIN            Plot with a minimum MIN [default: 16.0].
  --max=MAX -x MAX            Plot with a minimum MAX [default: 23.5].
  --X                         Plot with X set to the index. Default.
  --Y                         Plot with Y set to the index.
  --Z                         Plot with Z set to the index.
  --index=INDEX -i INDEX      The index that defines the cutplane. [default: 0].
  --half                      Just cut it in half along the specified axis.
  --T                         Transpose the SP.
  --no-log                    Do not log10 the data.
  --title=TITLE               Set title.
'''

import numpy as np;
import matplotlib;
import cPickle as pickle;
import math;
from docopt import docopt;
from misc import read;

def main():
    opts=docopt(__doc__,help=True);
    vmin = float(opts['--min']);
    vmax = float(opts['--max']);
    if not opts['--X'] and not opts['--Y'] and not opts['--Z']:
        opts['--X']=True;
    outfile = opts['<outfile>'];
    
    if outfile: matplotlib.use('Agg');
    import matplotlib.pyplot as plt;
        
    f = read(opts['<infile>'],dumpfull=True);
    S = f['s'];
    coords = [i for i in f.keys() if i in ['x','y','z']];
    if coords == ['x','y','z']:
        SP,label=prep3d(S);
        coords.remove(label);
        xlabel = f[coords[0]];
        ylabel = f[coords[1]];
        x = f[coords[0]];
        y = f[coords[1]];
    else:
        if len(coords) < 2:
            raise RuntimeError("Cannot make a 2D plot that isn't 2D data.")
        xlabel = f['0th'];
        ylabel = f['1st'];
        SP = S;
        x = f[xlabel];
        y = f[ylabel];
    #other preparation
    SP = np.nan_to_num(SP);
    if not opts['--no-log']: SP = np.log10(SP+0.1);
    if opts['--T']:
        SP=SP.T;
        t=x; x=y; y=t;
        t=xlabel; xlabel=ylabel; ylabel=t;
    xmin,xmax = x; ymin,ymax = y;
    #selecting index.                         v--This is intentional, to match matlab semantics.
    Y,X = np.mgrid[ ymin:ymax:len(SP[:,0])*1j,xmin:xmax:len(SP[0,:])*1j];
    plt.pcolormesh(X, Y, SP,vmin=vmin,vmax=vmax);
    plt.xlabel('{} ($\mu m$)'.format(xlabel));
    plt.ylabel('{} ($\mu m$)'.format(ylabel));
    c=plt.colorbar();
    c.set_label('log10 of density');
    if opts['--title']:
        plt.title(opts['--title']);
    if outfile:
        plt.savefig(opts['<outfile>']);
    else:
        plt.show();
pass;

def prep3d(S,opts):
    #selecting index.
    if opts['--index']:
        i = int(opts['--index']);
    elif opts['--half']:
        if opts['--X']:
            i = len(S[:,0,0])/2;
            label = 'x';
        elif opts['--Y']:
            i = len(S[0,:,0])/2;
            label = 'y';
        else:
            i = len(S[0,0,:])/2;
            label = 'z';
    else:
        i = 0;
    #selecting plane.
    if opts['--X']:
        SP = S[i,:,:];
    elif opts['--Y']:
        SP = S[:,i,:];
    else:
        SP = S[:,:,i];
    if opts['--T']: SP=SP.T
    return SP, label;

if __name__ == '__main__':
    main();
