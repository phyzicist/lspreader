#!/usr/bin/env python2
'''
Plot the angular/energy/charge plot.

Usage:
  angular.py [options] <input> <output>
  angular.py [options] <input>

Options:
  --polar-bins=BINS -p BINS   Set the number of polar bins.
  --radial-bins=BINS -r BINS  Set the number of radial bins.
  --title=TITLE -t Title      Set the title.
  --clabel=CLABEL             Set colorbar label.
  --KeV                       Scale by 100s of KeV instead of MeV.
  --max-E=MAXE                Set the maximum energy value (in units depending on --KeV flag
  --E-step=ESTEP              Set the step of grid lines for Energy.
'''
import numpy as np;
import matplotlib.pyplot as plt;
import cPickle
from matplotlib import colors;
from docopt import docopt;
import math as m;
from misc import conv;
def mk_cmap():
    hsv = np.array([[[0.00, 0.6, 1.0],
                     [0.10, 0.6, 1.0],
                     [0.20, 0.6, 1.0],
                     [0.30, 0.6, 1.0],
                     [0.40, 0.6, 1.0],
                     [0.50, 0.6, 1.0],
                     [0.60, 0.6, 1.0],
                     [0.70, 0.6, 1.0],
                     [0.80, 0.6, 1.0]]]);
    rgb = colors.hsv_to_rgb(hsv);
    def mk_component(cmp):
        l = len(cmp);
        inter = np.linspace(0.001,1.0,l);
        ret = [[i,j,j]  for i,j in zip(inter,cmp)];
        ret[0][1]=1.0;
        ret = [[0.0,1.0,1.0]]+ret;
        return tuple(ret);
    r = np.array(mk_component(rgb.T[0]));
    g = np.array(mk_component(rgb.T[1]));
    b = np.array(mk_component(rgb.T[2]));
    cd={'red':r,'green':g,'blue':b};
    return colors.LinearSegmentedColormap('cmap',cd, 1024);

def main():
    opts = docopt(__doc__,help=True);
    inname = opts['<input>'];
    outname = opts['<output>'];
    phi_spacing = conv(opts['--polar-bins'],default=180,func=float);
    r_spacing = conv(opts['--radial-bins'],default=40,func=float);
    clabel = conv(opts['--clabel'],default='$p C$');
    maxE = conv(opts['--max-E'],default=(400 if opts['--KeV'] else 4.0),func=float);
    rstep = conv(opts['--E-step'],default=(100 if opts['--KeV'] else 1.0),func=float)
    with open(inname,'r') as f:
        d = cPickle.load(f)
    r = d['KE'];
    if opts['--KeV']:
        r/=1e3;
    else:
        r/=1e6;
    s =-d['q']*1e6;
    plt.hist(s,bins=50);
    phi_bins = np.linspace(min(d['phi']),max(d['phi']),phi_spacing+1);
    r_bins   = np.linspace(0, maxE, r_spacing+1);
    PHI,R = np.mgrid[ min(d['phi']) : max(d['phi']) : phi_spacing*1j,
                      0 : maxE : r_spacing*1j];
    S,_,_ = np.histogram2d(d['phi'],r,bins=(phi_bins,r_bins),weights=s);
    fig = plt.figure(1);
    ax=plt.subplot(projection='polar',axisbg='white');
    #surf=ax.pcolormesh(PHI,R,S,vmin=s.min(),vmax=s.max(),cmap=mk_cmap());
    surf=plt.pcolormesh(PHI,R,S,cmap=mk_cmap());
    #making radial guides. rgrids only works for plt.polar calls
    full_phi = np.linspace(0.0,2*np.pi,100);
    for i in np.arange(0.0,maxE,rstep)[1:]:
        plt.plot(full_phi,np.ones(full_phi.shape)*i,c='gray', lw=1, ls='--');
    ax.set_theta_zero_location('S');
    unit = 'KeV' if opts['--KeV'] else 'MeV';
    label_str = '{} '+unit;
    labels    = np.arange(0.0,maxE,rstep)[1:];
    plt.rgrids(labels,labels=map(label_str.format,labels),angle=190);
    ax.set_xticklabels(['Backwards\n0$^{\circ}$',
                        '45$^{\circ}$',
                        'left\n90$^{\circ}$',
                        '135$^{\circ}$',
                        'Forwards\n180$^{\circ}$',
                        '215$^{\circ}$',
                        'right\n270$^{\circ}$',
                        '315$^{\circ}$']);
    c=fig.colorbar(surf);
    c.set_label(clabel);
    if opts['--title']:
        plt.title(opts['--title'],loc='left',fontdict={'fontsize':28});
    if outname:
        plt.savefig(outname);
    else:
        plt.show();
    pass;

if __name__ == "__main__":
    main();
