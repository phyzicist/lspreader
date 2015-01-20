#!/usr/bin/env python2
'''
Plot the angular/energy/charge plot.

Usage:
  angular.py [options] <input> <output>
  angular.py [options] <input>

Options:
  --angle-bins=BINS -a BINS   Set the number of angle bins.  [default: 180]
  --radial-bins=BINS -r BINS  Set the number of radial bins. [default: 40]
  --title=TITLE -t Title      Set the title. 
  --clabel=CLABEL -c CLABEL   Set colorbar label. [default: $p C$]
  --KeV -k                    Scale by 100s of KeV instead of MeV.
  --max-e=MAXE -e MAXE        Set the maximum energy value (in units depending on --KeV flag)
  --e-step=ESTEP              Set the step of grid lines for Energy.
  --high-res -H               Output a high resolution plt.
  --max-q=MAXQ -q MAXQ        Set the maximum for the charge (pcolormesh's vmax value).
  --normalize -n              Normalize the histogram to 1 *eV^-1 rad^-1 .
  --factor=F -f F             Multiply histogram by F. [default: 1.0]
  --polar -p                 Plot polar angles, letting the east direction be forward.
'''
import numpy as np;
import matplotlib.pyplot as plt;
import cPickle as pickle;
from matplotlib import colors;
from docopt import docopt;
from misc import conv,pastel;

def main():
    opts = docopt(__doc__,help=True);
    inname = opts['<input>'];
    outname = opts['<output>'];
    phi_spacing = float(opts['--angle-bins']);
    E_spacing = float(opts['--radial-bins']);
    clabel = opts['--clabel'];
    maxE = conv(opts['--max-e'],default=(1000 if opts['--KeV'] else 4.0),func=float);
    maxQ = float(opts['--max-q']) if opts['--max-q'] else None;
    Estep = conv(opts['--e-step'],default=(250 if opts['--KeV'] else 1.0),func=float);
    F = float(opts['--factor']);
    with open(inname,'r') as f:
        d = pickle.load(f)
    e = d['KE'];
    if opts['--KeV']:
        e/=1e3;
    else:
        e/=1e6;
    s =-d['q']*1e6*F;
    if opts['--polar']:
        phi = d['phi_n'];
        phi_labels = ['Forward\n0$^{\circ}$',
                      '45$^{\circ}$',
                      'Up\n90$^{\circ}$',
                      '135$^{\circ}$',
                      'Backwards\n180$^{\circ}$',
                      '215$^{\circ}$',
                      'Down\n270$^{\circ}$',
                      '315$^{\circ}$'];
    else:
        phi = d['phi'];
        phi_labels = ['Forward\n0$^{\circ}$',
                      '45$^{\circ}$',
                      'Left\n90$^{\circ}$',
                      '135$^{\circ}$',
                      'Backwards\n180$^{\circ}$',
                      '215$^{\circ}$',
                      'Right\n270$^{\circ}$',
                      '315$^{\circ}$'];
    phi_bins = np.linspace(-np.pi,np.pi,phi_spacing+1);
    E_bins   = np.linspace(0, maxE, E_spacing+1);
    PHI,E = np.mgrid[ -np.pi : np.pi : phi_spacing*1j,
                      0 : maxE : E_spacing*1j];
    S,_,_ = np.histogram2d(phi,e,bins=(phi_bins,E_bins),weights=s);
    if opts['--normalize']:
        Efactor = maxE/E_spacing;
        if opts['--KeV']:
            Efactor *= 1e-3;
            clabel += ' rad$^{-1}$ KeV$^{-1}$'
        else:
            clabel += ' rad$^{-1}$ MeV$^{-1}$'
        S /= Efactor * 2*np.pi/phi_spacing;
    fig = plt.figure(1);
    ax=plt.subplot(projection='polar',axisbg='white');
    surf=plt.pcolormesh(PHI,E,S,cmap=pastel,vmax=maxQ);
    #making radial guides. rgrids only works for plt.polar calls
    full_phi = np.linspace(0.0,2*np.pi,100);
    for i in np.arange(0.0,maxE,Estep)[1:]:
        plt.plot(full_phi,np.ones(full_phi.shape)*i,c='gray', lw=1, ls='--');
    ax.set_theta_zero_location('N');
    unit = 'KeV' if opts['--KeV'] else 'MeV';
    rlabel_str = '{} '+unit;
    rlabels    = np.arange(0.0,maxE,Estep)[1:];
    plt.rgrids(rlabels, labels=map(rlabel_str.format,rlabels),angle=350);
    ax.set_xticklabels(phi_labels);
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
