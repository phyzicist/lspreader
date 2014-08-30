#!/usr/bin/env python2
'''
Choose particles positions for tracks.

Usage:
  pmovie_select.py [options] <input> <output> <select>
  pmovie_select.py [options] --test-only <input>

Options:
  --pool-size PS -p  PS      Set the size of the processor pool.
  --xilim=XILIM              Set initial X limits, using a python list as a string.
  --yilim=YILIM              Set initial Y limits, using a python list as a string.
  --zilim=ZILIM              Set initial Z limits, using a python list as a string.
  --select=N                 Randomly select N initial particle positions.
'''
from docopt import docopt;
import lspreader as rd;
import numpy as np;
from misc import conv;

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
        print(f.header);
        data=f._getmovie(ps);
    if opts['--count-only']:
        print(len(d[0]['x']));
    else:
        pass;
        #here we play a little :) We randomly choose indices, and see if they fall inside.
        #there are some assumptions here:
        #that there are only a few ipp's, thi
if __name__ == '__main__':
    main();
