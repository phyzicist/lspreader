#!/usr/bin/env python2
'''
Make pmovie particles.

Usage:
  pmovie.py [options] <select-file> <input>

Options:
  --parallel=PS -p PS       Set parallel with pool size. [default: 1]
  --subdivisions=S -s S     Run sequentially with the following subdivisions S of the pmovie. [default: 8]
  --deep-compare -d         Run with a more correct deep checking.
  --verbose -v              Run verbosely.
'''
import lspreader as rd;
import cPickle;
import numpy as np;
from docopt import docopt;
import random;
import itertools as it;
from time import time;
import multiprocessing;
import sys;

def _print(s):
    print(s);
    sys.stdout.flush();

st=time();
opts = docopt(__doc__,help=True);
inname = opts['<input>'];
selectfile = opts['<select-file>'];
ss = int(opts['--subdivisions'])
ps = int(opts['--parallel']);
vprint = _print if opts['--verbose'] else lambda s: None;

vprint('reading in select file {}'.format(selectfile));

with open(selectfile,'r') as f:
    s = cPickle.load(f);
selecti = np.array([s['xi'],s['yi'],s['zi']]).T;
Ns = len(selecti);
xlim = s['xlim'];
ylim = s['ylim'];
zlim = s['zlim'];
del s;

vprint('reading in {}'.format(inname));
with rd.LspOutput(inname,verbose=opts['--verbose']) as f:
    d=f._getmovie();

frame = d[0];
keys = [k for k in frame if k not in ['t','step','pnum']];
del d;
d=frame['data'];
del frame['data'];
vprint('filtering');
goods = (xlim[0] <= d['xi']) & (d['xi'] <= xlim[1]) & \
        (ylim[0] <= d['yi']) & (d['yi'] <= ylim[1]) & \
        (zlim[0] <= d['zi']) & (d['zi'] <= zlim[1]);

d = d[goods];
del goods

if opts['--deep-compare']:
    match = lambda x,y,z: np.isclose(x,d['xi']) & np.isclose(y,d['yi']) & np.isclose(z,d['zi']);
else:
    #this assumes the ipp's are copied, not calculated.
    match = lambda x,y,z: (x==d['xi']) & (y==d['yi']) & (z==d['zi']);

def getmatch(i):
    return np.argmax(match(i[0],i[1],i[2]))

vprint('starting reading');
pool = multiprocessing.Pool(ps);
results = pool.map(getmatch,selecti);
pool.close();
results = np.array(results);
d = d[results];
frame['data'] = d;
frame['valid'] = results != 0;
frame['valid'][0] = True;

outname = "points{}.pt".format(frame['step']);
vprint('outputting {}'.format(outname));
with open(outname,'w') as f:
    cPickle.dump(frame,f,2);
pass;
vprint('finished with time {} min'.format((time()-st)/60.0));



