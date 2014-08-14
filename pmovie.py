#!/usr/bin/env python2
'''
Choose particles positions for tracks.

Usage:
  pmovie.py [options] <output> <input>...

Options:
  --xilim=XILIM         Set initial X limits, using a python list as a string.
  --yilim=YILIM         Set initial Y limits, using a python list as a string.
  --zilim=ZILIM         Set initial Z limits, using a python list as a string.
  --skip=SKIP           Set skipping interval.
  --pool-size=PS        Set pool size.
  --select=N            Randomly select N trajectories.
'''
import lspreader as rd;
import cPickle;
import numpy as np;
from docopt import docopt;
import re;
import random;

def main():
    opts = docopt(__doc__,help=True);
    outname = opts['<output>'];
    names = opts['<input>'];
    ps = int(opts['--pool-size']) if opts['--pool-size'] else 8;
    skip = int(opts['--skip']) if opts['--skip'] else 1;
    inflim = (float('-inf'),float('inf'));
    #Pretty unsafe. This is a script, not a cryptographically secure app.
    xlim = eval(opts['--xilim']) if opts['--xilim'] else inflim;
    ylim = eval(opts['--yilim']) if opts['--yilim'] else inflim;
    zlim = eval(opts['--zilim']) if opts['--zilim'] else inflim;
    sel  = int(opts['--select']) if opts['--select'] else 100;
    frames=[];
    for name in names:
        print('reading in {}'.format(name));
        with rd.LspOutput(name) as f:
            frames.extend(f._getmovie(ps,skip));
    print('setting up first frame');
    data={};
    for p in frames[0]['particles']:
        pos = (p['xi'],p['zi']);
        if  xlim[0] <= pos[0] <= xlim[1] and \
            zlim[0] <= pos[1] <= zlim[1]:#ylim[0] <= pos[1] <= ylim[1] and 
            pos=str(pos);
            data[pos]=[[p['x']],[p['z']]];
    print('reading frames');
    for frame in frames[1:]:
        for p in frame['particles']:
            pos = (p['xi'],p['zi']);
            if  xlim[0] <= pos[0] <= xlim[1] and \
                zlim[0] <= pos[1] <= zlim[1]:#ylim[0] <= pos[1] <= ylim[1] and \
                pos=str(pos);
                data[pos][0].append(p['x']);
                #data[pos][1].append(p['y']);
                data[pos][1].append(p['z']);
        del frame['particles'];
    pass;
    del frames;
    #finally, selecting
    print('selecting {} trajectories'.format(sel));
    out = random.sample([data[k] for k in data],sel);
    del data;
    print('outputting {}'.format(outname));
    with open(outname,'w') as f:
        cPickle.dump(out,f,2);
    pass;

if __name__ == '__main__':
    main();
