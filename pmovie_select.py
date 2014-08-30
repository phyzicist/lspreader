#!/usr/bin/env python2
'''
Choose particles positions for tracks.

Usage:
  pmovie_select.py [options] <input> <output>
  pmovie_select.py [options] --count-only <input>

Options:
  --pool-size PS -p  PS      Set the size of the processor pool.
  --xilim=XILIM              Set initial X limits, using a python list as a string.
  --yilim=YILIM              Set initial Y limits, using a python list as a string.
  --zilim=ZILIM              Set initial Z limits, using a python list as a string.
  --select=N                 Randomly select N initial particle positions [default: 100].
'''
from docopt import docopt;
import lspreader as rd;
import numpy as np;
import random as rnd;
from misc import conv;
import cPickle;

def main():
    opts = docopt(__doc__,help=True);
    outname = opts['<output>'];
    name = opts['<input>'];
    ps = conv(opts['--pool-size'],default=8,func=float);
    inflim = (float('-inf'),float('inf'));
    xlim = conv(opts['--xilim'],default=inflim,func=eval);
    ylim = conv(opts['--yilim'],default=inflim,func=eval);
    zlim = conv(opts['--zilim'],default=inflim,func=eval);
    select = conv(opts['--select'],default=100,func=int);
    print('reading in {}'.format(name));
    with rd.LspOutput(name,verbose=True) as f:
        d=f._getmovie(ps);
    frame = d[0];
    del d[:1];
    N = len(frame['x']);
    if opts['--count-only']:
        print('There are {} particles.'.format(N));
        exit(0);
    #here we play a little :) We randomly choose indices, and see if they fall inside.
    #The main assumption here is that there are only a few ipp that are being selected.
    #that there are only a few ipp's that we are selecting
    out = {};
    keys = [k for k in frame if k not in ['t','step']];
    for k in keys:
        out[k] = [];
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
        pass;
    print('found {} points out of {} trials'.format(len(out['x']),
                                                    len(ns)));
    print('outputting to {}'.format(outname));
    with open(outname,'wb') as f:
        cPickle.dump(out,f,2);
    pass;

if __name__ == '__main__':
    main();
