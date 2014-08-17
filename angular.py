#!/usr/bin/env python2
'''
Plot the angular/energy/charge plot.

Usage:
  angular.py [options] <input> <output>
  angular.py [options] <input>

Options:
  --polar-bins=BINS -p BINS   Set the number of polar bins.
  --radial-bins=BINS -r BINS  Set the number of radial bins.
'''
import numpy as np;
import matplotlib.pyplot as plt;
import cPickle
from matplotlib import colors;
from docopt import docopt;

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
    phi_spacing = float(opts['--polar-bins']) if opts['--polar-bins'] else 180;
    r_spacing   = float(opts['--radial-bins']) if opts['--radial-bins'] else 40;
    
    with open(inname,'r') as f:
        d = cPickle.load(f)
    r = d['KE']/1e6;
    s =-d['q']*1e6;
    plt.hist(s,bins=50);
    phi_bins = np.linspace(min(d['phi']),max(d['phi']),phi_spacing+1);
    r_bins   = np.linspace(0, 4.0, r_spacing+1);
    PHI,R = np.mgrid[ min(d['phi']) : max(d['phi']) : phi_spacing*1j,
                      0 : 4.0 : r_spacing*1j];
    S,_,_ = np.histogram2d(d['phi'],r,bins=(phi_bins,r_bins),weights=s);
    fig = plt.figure(1);
    ax = plt.subplot(projection='polar',axisbg='white');
    #surf=ax.pcolormesh(PHI,R,S,vmin=s.min(),vmax=s.max(),cmap=mk_cmap());
    surf=ax.pcolormesh(PHI,R,S,cmap=mk_cmap());
    ax.plot(PHI,np.ones(PHI.shape)*1,c='gray', lw=1, ls='--');
    ax.plot(PHI,np.ones(PHI.shape)*2,c='gray', lw=1, ls='--');
    ax.plot(PHI,np.ones(PHI.shape)*3,c='gray', lw=1, ls='--');
    ax.set_theta_zero_location('S');
    ax.set_yticks(np.arange(0.0,4.0));
    ax.annotate('Energy (MeV)',xy=(np.pi/8+0.1,1.5),rotation=-67.5);
    c=fig.colorbar(surf);
    c.set_label('$p C$');
    if outname:
        plt.savefig(outname);
    else:
        plt.show();
    pass;

if __name__ == "__main__":
    main();
