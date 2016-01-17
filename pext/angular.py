#!/usr/bin/env python2
'''
Show an angular/energy/charge plot movie.

Usage:
  angular.py [options] <input>
  angular.py [options] <input> <output>

Options:
  --angle-bins=BINS -a BINS   Set the number of angle bins.  [default: 180]
  --radial-bins=BINS -r BINS  Set the number of radial bins. [default: 40]
  --ltitle=TITLE -t TITLE     Set the title.
  --rtitle=TITLE -T TITLE     Set the right title.
  --clabel=CLABEL -c CLABEL   Set colorbar label. [default: $p C$]
  --no-cbar                   Turn off the colorbar.
  --KeV -k                    Scale by 100s of KeV instead of MeV.
  --max-e=MAXE -e MAXE        Set the maximum energy value (in units depending on --KeV flag)
  --e-step=ESTEP              Set the step of grid lines for Energy.
  --high-res -H               Output a high resolution plt.
  --max-q=MAXQ -q MAXQ        Set the maximum for the charge (pcolormesh's vmax value).
  --min-q=MINQ                Set a minimum charge.
  --normalize -n              Normalize the histogram to 1 *eV^-1 rad^-1 .
  --factor=F -f F             Multiply histogram by F. [default: 1.0]
  --polar -p                  Plot polar angles, letting the east direction be forward.
  --oap=ANGLE -o ANGLE        Set the width angle of the OAP. [default: 50.47]
  --log10 -l                  Plot a logarithmic pcolor instead of a linear one.
  --cmap=CMAP                 Use the following cmap [default: pastel].
  --e-direction=ANGLE         The angle for the radial labels.
  --e-units=UNIT              The units for the radial labels.
  --agg                       Use the agg backend.
'''
import numpy as np;
import matplotlib;
import matplotlib.pyplot as plt;
import matplotlib.patheffects as pe;
import cPickle as pickle;
from matplotlib import colors;
from docopt import docopt;
from misc import conv, test, readfile
from cmaps import pastel_clear,plasma_clear,viridis_clear,magma_clear_r;
import re;

def prep(opts):
    '''I put this here in order to reuse this'''
    inname = opts['<input>'];
    outname = opts['<output>'];
    
    phi_spacing = float(opts['--angle-bins']);
    E_spacing = float(opts['--radial-bins']);    
    maxE = conv(opts['--max-e'],default=(1000 if opts['--KeV'] else 4.0),func=float);
    maxQ = float(opts['--max-q']) if opts['--max-q'] else None;
    Estep = conv(opts['--e-step'],default=(250 if opts['--KeV'] else 1.0),func=float); 
    F = float(opts['--factor']);
    d = np.load(inname, allow_pickle=True);
    
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
    kw = {
        'angle_bins' : float(opts['--angle-bins']),
        'radial_bins': float(opts['--radial-bins']),
        'max_e': float(opts['--max-e']) if opts['--max-e'] else (1000 if opts['--KeV'] else 4.0),
        'max_q': float(opts['--max-q']) if opts['--max-q'] else None,
        'min_q': float(opts['--min-q']) if opts['--min-q'] else None,
        'KeV': opts['--KeV'],
        'clabel' : opts['--clabel'],
        'colorbar' : not opts['--no-cbar'],
        'e_step' : float(opts['--e-step']) if opts['--e-step'] else None,
        'labels':phi_labels,
        'rtitle':opts['--rtitle'],
        'ltitle':opts['--ltitle'],
        'outname':outname,
        'oap': float(opts['--oap']) if opts['--oap'] != 'none' else None,
        'log_q': opts['--log10'],
    };
    if opts['--cmap'] == 'viridis':
        kw['cmap'] = viridis_clear;
    elif opts['--cmap'] == 'plasma':
        kw['cmap'] = plasma_clear;
    elif opts['--cmap'] == 'magma_r':
        kw['cmap'] = magma_clear_r;
    else:
        kw['cmap'] = pastel_clear;
    if opts['--normalize']:
        Efactor = kw['max_e']/kw['radial_bins'];
        if kw['KeV']:
            Efactor *= 1e-3;
            kw['clabel'] += ' rad$^{-1}$ MeV$^{-1}$'            
        else:
            kw['clabel'] += ' rad$^{-1}$ MeV$^{-1}$'
        s /= Efactor*2*np.pi/phi_spacing;
    kw['rgridopts'] = {};
    if opts['--e-direction']:
        kw['rgridopts'].update({'angle':opts['--e-direction']});
    if opts['--e-units']:
        kw['rgridopts'].update({'unit':opts['--e-units']});
    
    return s,phi,e,kw,d;

def main():
    opts=docopt(__doc__,help=True);
    s,phi,e,kw,_ = prep(opts);
    if opts['<output>'] and opts['--agg']:
        plt.change_backend('agg');
    angular(s,phi,e,**kw);
    if opts['<output>']:
        if opts['--high-res']:
            plt.savefig(opts['<output>'],dpi=1000);
        else:
            plt.savefig(opts['<output>']);
    else:
        plt.show();
    
    pass;


def angular(s, phi, e,
            colorbar=True,**kw):
    '''
    Make the angular plot.

    Arguments:
      s   -- the charges.
      phi -- the angles of ejection.
      e   -- energies of each charge.

    Keyword Arugments:
      max_e       -- Maximum energy.
      max_q       -- Maximum charge.
      angle_bins  -- Set the number of angle bins.
      radial_bins -- Set the number of radial (energy) bins.
      clabel      -- Set the colorbar label.
      colorbar    -- If true, plot the colorbar.
      e_step      -- Set the steps of the radius contours.
      labels      -- Set the angular labels.
      KeV         -- Use KeV isntead of MeV.
      fig         -- If set, use this figure, Otherwise,
                     make a new figure.
      ax          -- If set, use this axis. Otherwise,
                     make a new axis.
      ltitle      -- Make a plot on the top left.
      rtitle      -- Make a plot on the top right.
      log_q       -- log10 the charges.
      cmap        -- use the colormap cmap.
      rgridopts   -- pass a dictionary that sets details for the
                     rgrid labels.
    '''
    phi_spacing = kw['angle_bins'];
    E_spacing =   kw['radial_bins'];    
    maxE  = kw['max_e']  if kw['max_e'] else (1000 if kw['KeV'] else 4.0);
    maxQ  = kw['max_q']  if kw['max_q'] else None;
    minQ  = kw['min_q']  if kw['min_q'] else None;
    Estep = kw['e_step'] if kw['e_step'] else (250 if kw['KeV'] else 1.0);
    clabel = kw['clabel'] if kw['clabel'] else '$pC$';
    if test(kw,'cmap'):
        cmap = kw['cmap'];
    else:
        cmap = pastel_clear;
    phi_bins = np.linspace(-np.pi,np.pi,phi_spacing+1);
    E_bins   = np.linspace(0, maxE, E_spacing+1);
            
    PHI,E = np.mgrid[ -np.pi : np.pi : phi_spacing*1j,
                      0 : maxE : E_spacing*1j];
    S,_,_ = np.histogram2d(phi,e,bins=(phi_bins,E_bins),weights=s);
    fig = kw['fig'] if test(kw,'fig') else plt.figure(1,facecolor=(1,1,1));
    ax  = kw['ax'] if test(kw,'ax') else plt.subplot(projection='polar',axisbg='white');
    norm = matplotlib.colors.LogNorm() if test(kw,'log_q') else None;
    
    surf=plt.pcolormesh(PHI,E,S,norm=norm, cmap=cmap,vmin=minQ,vmax=maxQ);
    #making radial guides. rgrids only works for plt.polar calls
    full_phi = np.linspace(0.0,2*np.pi,100);
    for i in np.arange(0.0,maxE,Estep)[1:]:
        plt.plot(full_phi,np.ones(full_phi.shape)*i,
                 c='gray',alpha=0.9,
                 lw=1, ls='--');
    ax.set_theta_zero_location('N');
    ax.patch.set_alpha(0.0);
    ax.set_axis_bgcolor('red');
    #making rgrid
    if test(kw, 'rgridopts'):
        ropts = kw['rgridopts'];
        if test(ropts, 'unit'):
            runit = ropts['unit'];
        else:
            runit = 'KeV' if test(kw,'KeV') else 'MeV';
        if test(ropts, 'angle'):
            rangle = ropts['angle'];
        else:
            rangle = 45;
        if test(ropts, 'size'):
            rsize = ropts['size'];
        else:
            rsize = 10.5;
        if test(ropts, 'invert'):
            c1,c2 = "w","black";
        else:
            c1,c2 = "black","w";
    else:
        runit = 'KeV' if test(kw,'KeV') else 'MeV';
        rangle = 45;
        rsize = 10.5;
        c1,c2 = "black","w";
    rlabel_str = '{} ' + runit;
    rlabels    = np.arange(0.0,maxE,Estep)[1:];
    _,ts=plt.rgrids(rlabels,
                    labels=map(rlabel_str.format,rlabels),
                    angle=rangle);
    for t in ts:
        t.set_path_effects([
            pe.Stroke(linewidth=1.5, foreground=c2),
            pe.Normal()
        ]);
        t.set_size(rsize);
        t.set_color(c1);
    if test(kw,'oap'):
        oap = kw['oap']/2 * np.pi/180;
        maxt = oap+np.pi; mint = np.pi-oap;
        maxr  = maxE*.99;
        if kw['KeV']:
            minr=120;
        else:
            minr=.12;
        ths=np.linspace(mint, maxt, 20);
        rs =np.linspace(minr, maxr, 20);
        mkline = lambda a,b: plt.plot(a,b,c=(0.2,0.2,0.2),ls='-',alpha=0.5);
        mkline(ths, np.ones(ths.shape)*minr)
        mkline(mint*np.ones(ths.shape), rs);
        mkline(maxt*np.ones(ths.shape), rs);
    if test(kw,'labels'):
        ax.set_xticks(np.pi/180*np.linspace(0,360,len(kw['labels']),endpoint=False));
        ax.set_xticklabels(kw['labels']);
    if colorbar:
        c=fig.colorbar(surf,pad=0.1);
        c.set_label(clabel);
    if test(kw,'ltitle'):
        if len(kw['ltitle']) <= 4:
            ax.set_title(kw['ltitle'],loc='left',fontdict={'fontsize':28});
        else:
            ax.text(np.pi/4+0.145,maxE+Estep*2.5,kw['ltitle'],fontdict={'fontsize':28});
    if test(kw,'rtitle'):
        if '\n' in kw['rtitle']:
            fig.text(0.60,0.875,kw['rtitle'],fontdict={'fontsize':22});
        else:
            plt.title(kw['rtitle'],loc='right',fontdict={'fontsize':22});
    return (surf, ax, fig, (phi_bins, E_bins));

if __name__ == "__main__":
    main();
