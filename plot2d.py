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
  --T -T                      Transpose the SP.
  --no-log                    Do not log10 the data.
  --title=TITLE -t TITLE      Set title.
  --highlight=VAL             Highlight a value on the cmap.
  --high-res -H               Use high resolution.
  --clabel=CLABEL -c CLABEL   Set the colorbar label. [default: log10 of electron density (cm$^{-3}$)].
  --fontsize SIZE -f SIZE     Set the fontsize. [default: 14].
'''

import numpy as np;
import matplotlib;
import cPickle as pickle;
import math;
from docopt import docopt;
import misc as m;

def main():
    opts=docopt(__doc__,help=True);
    vmin = float(opts['--min']);
    vmax = float(opts['--max']);
    if not opts['--X'] and not opts['--Y'] and not opts['--Z']:
        opts['--X']=True;
    outfile = opts['<outfile>'];
    
    if outfile: matplotlib.use('Agg');
    import matplotlib.pyplot as plt;
        
    f = m.readfile(opts['<infile>'],dumpfull=True);
    
    S = f['s'];
    coords = [i for i in f.keys() if i in ['x','y','z']];
    coords.sort();
    if coords == ['x','y','z']:
        SP, label = prep3d(S, opts);
        coords.remove(label);
        xlabel = coords[1];
        ylabel = coords[0];
        x = f[coords[1]];
        y = f[coords[0]];
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
    #convert to microns
    xmin*=1e4;
    xmax*=1e4;
    ymin*=1e4;
    ymax*=1e4;
    #selecting index.                         v--This is intentional, to match matlab semantics.
    Y,X = np.mgrid[ ymin:ymax:len(SP[:,0])*1j,xmin:xmax:len(SP[0,:])*1j];
    if opts['--highlight']:
        val = float(opts['--highlight']);
        cmap = m.mkstrip_cmap(vmin,vmax,val,[0.6,0.3,0.0]);
    else:
        cmap = m.pastel;
    plt.pcolormesh(X, Y, SP,vmin=vmin,vmax=vmax,cmap=cmap);
    plt.xlim(xmin,xmax);
    plt.ylim(ymin,ymax);
    fontsize = float(opts['--fontsize']);
    plt.xlabel('{} ($\mu m$)'.format(xlabel),{'fontsize':fontsize});
    plt.ylabel('{} ($\mu m$)'.format(ylabel),{'fontsize':fontsize});
    c=plt.colorbar();
    c.set_label(opts['--clabel'],size=fontsize);
    c.ax.tick_params(labelsize=fontsize-2);
    if opts['--title']:
        plt.title(opts['--title'],{'fontsize':fontsize+4});
    if outfile:
        if opts['--high-res']:
            plt.savefig(outfile,dpi=1000);
        else:
            plt.savefig(outfile);
    else:
        plt.show();
pass;

def prep3d(S,opts):
    #selecting index.
    if opts['--half']:
        if opts['--X']:
            i = len(S[:,0,0])/2;
        elif opts['--Y']:
            i = len(S[0,:,0])/2;
        else:
            i = len(S[0,0,:])/2;
    elif opts['--index']:
        i = int(opts['--index']);
    else:
        i = 0;
    #selecting plane.
    if opts['--X']:
        SP = S[i,:,:];
        label = 'x';
    elif opts['--Y']:
        SP = S[:,i,:];
        label = 'y';
    else:
        SP = S[:,:,i];
        label = 'z';
    if opts['--T']: SP=SP.T
    return SP, label;

if __name__ == '__main__':
    main();
