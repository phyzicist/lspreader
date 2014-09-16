#!/usr/bin/env python2
'''
Choose particles positions for tracks.

Usage:
  pmovie_select.py [options] <input> <output>
  pmovie_select.py [options] --count-only <input>
  pmovie_select.py [options] --histograms <input> <output>

Options:
  --pool-size PS -p  PS      Set the size of the processor pool [default: 8].
  --xilim=XILIM              Set initial X limits, using a python list as a string.
  --yilim=YILIM              Set initial Y limits, using a python list as a string.
  --zilim=ZILIM              Set initial Z limits, using a python list as a string.
  --select=N                 Randomly select N initial particle positions [default: 100].
  --verbose -v               Allow verbose output.
'''
from docopt import docopt;
import lspreader as rd;
import numpy as np;
import random as rnd;
from misc import conv;
import cPickle;

def logprint(s):
    print(s);

def main():
    opts = docopt(__doc__,help=True);
    outname = opts['<output>'];
    name = opts['<input>'];
    ps = conv(opts['--pool-size'],default=8,func=int);
    inflim = (float('-inf'),float('inf'));
    xlim = conv(opts['--xilim'],default=inflim,func=eval);
    ylim = conv(opts['--yilim'],default=inflim,func=eval);
    zlim = conv(opts['--zilim'],default=inflim,func=eval);
    select = conv(opts['--select'],default=100,func=int);
    if opts['--verbose']:
        vprint = logprint;
    else:
        vprint = lambda s: None;
    vprint('reading in {}'.format(name));
    with rd.LspOutput(name,verbose=True) as f:
        d=f._getmovie(ps);
    frame = d[0];
    del d[:1];
    N = frame['pnum'];
    if opts['--count-only']:
        print('There are {} particles.'.format(N));
        exit(0);
    elif opts['--histograms']:
        import matplotlib;
        matplotlib.use('Agg');
        import matplotlib.pyplot as plt;
        plt.subplot(311);
        plt.title('xi');
        plt.hist(frame['xi'],bins=100);
        plt.subplot(312);
        plt.title('yi');
        plt.hist(frame['yi'],bins=100);
        plt.subplot(313);
        plt.title('zi');
        plt.hist(frame['zi'],bins=100);
        plt.savefig(outname);
        exit(0);
    #here we play a little :) We randomly choose indices, and see if they fall inside.
    #The main assumption here is that there are only a few ipp that are being selected.
    #that there are only a few ipp's that we are selecting
    keys = [k for k in frame if k not in ['t','step','pnum']];
    out = { k:[] for k in keys}
    ns = [];
    while len(out['x']) < select:
        i=rnd.randint(0,N-1);
        if i in ns:
            continue;
        ns.append(i);
        if  xlim[0] <= frame['xi'][i] <= xlim[1] and \
            ylim[0] <= frame['yi'][i] <= ylim[1] and \
            zlim[0] <= frame['zi'][i] <= zlim[1]:
            for k in out:
                out[k].append(frame[k][i]);
            vprint('now, we have {} particles'.format(len(out['x'])));
        pass;
    vprint('found {} points out of {} trials'.format(len(out['x']),
                                                     len(ns)));
    #making indices
    out['i'] = range(len(out['x']));
    vprint('outputting to {}'.format(outname));
    with open(outname,'wb') as f:
        cPickle.dump(out,f,2);
    pass;

if __name__ == '__main__':
    main();
