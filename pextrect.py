#!/usr/bin/env python2
'''
Plot a charge vs. two other components, x and y of a pext.py output.

Usage:
  angular.py [options] <input> <x> <y> [<output>]

Options:
  --x-bins=BINS -x BINS       Set the number of x bins. [default: 100]
  --y-bins=BINS -y BINS       Set the number of y bins. [default: 100]
  --x-label=LABEL LABEL       Set the x label.
  --y-label=LABEL LABEL       Set the y label.
  --xlim=LIM                  Set limits on the x variable.
  --ylim=LIM                  Set limits on the y variable.
  --x-no-restrict             Set limits, but don't restrict x.
  --y-no-restrict             Set limits, but don't restrict y.
  --title=TITLE -t TITLE      Set the title.
  --clabel=CLABEL             Set colorbar label. [default: $p C$]
  --high-res -H               Output a high resolution plt.
  --max-Q=MAXQ                Set the maximum for the charge (pcolormesh's vmax value).
  --normalize                 Normalize the histogram to 1 *eV^-1 rad^-1 .
  --factor=F                  Multiply histogram by F. [default: 1.0]
'''
import numpy as np;
import matplotlib.pyplot as plt;
import cPickle as pickle;
from docopt import docopt;
from misc import conv,pastel_rainbow;

def restrict(x,xlim):
    good = (x >= xlim[0]) & (x <= xlim[1]);
    return x[good];

def main():
    opts = docopt(__doc__,help=True);
    inname = opts['<input>'];
    outname = opts['<output>'];

    xname = opts['<x>'];
    yname = opts['<y>'];
    
    x_spacing = float(opts['--x-bins']);
    y_spacing = float(opts['--y-bins']);

    x_label = conv(opts['--x-label'],xname);
    y_label = conv(opts['--y-label'],yname);
    
    clabel = opts['--clabel'] if opts['--clabel'] else '$p C$';
    maxQ = float(opts['--max-Q']) if opts['--max-Q'] else None;
    F = float(opts['--factor']);
    with open(inname,'r') as f:
        d = pickle.load(f)

    x = d[xname];
    y = d[yname];
    if opts['--xlim'] and not opts['--x-no-restrict']:
        xlim = eval(opts['--xlim']);
        x = restrict(x,xlim);
    else:
        xlim = (x.min(), x.max());

    if opts['--ylim'] and not opts['--y-no-restrict']:
        ylim = eval(opts['--ylim']);
        y = restrict(y,ylim);
    else:
        ylim = (y.min(), y.max());

    s =-d['q']*1e6*F;
    
    x_bins = np.linspace(xlim[0],xlim[1],x_spacing+1);
    y_bins = np.linspace(ylim[0],ylim[1],y_spacing+1);
    
    X,Y = np.mgrid[ xlim[0] : xlim[1] : x_spacing*1j,
                    ylim[0] : ylim[1] : y_spacing*1j];
    S,_,_ = np.histogram2d(x,y,bins=(x_bins,y_bins),weights=s);
    if opts['--normalize']:
        S /= np.abs(xlim[1]-xlim[0])/x_spacing;
        S /= np.abs(ylim[1]-ylim[0])/y_spacing;
    fig = plt.figure(1);
    plt.xticks(x_bins[::x_spacing/10]);
    plt.yticks(y_bins[::y_spacing/10]);
    plt.xlabel(x_label);
    plt.ylabel(y_label);
    surf=plt.pcolormesh(X,Y,S,cmap=pastel_rainbow,vmax=maxQ);
    c=fig.colorbar(surf,pad=0.075);
    c.set_label(clabel);
    if opts['--title']:
        plt.title(opts['--title'],loc='left',fontdict={'fontsize':28});
    if outname:
        if opts['--high-res']:
            plt.savefig(outname,dpi=1000);
        else:
            plt.savefig(outname);
    else:
        plt.show();
    pass;

if __name__ == "__main__":
    main();
