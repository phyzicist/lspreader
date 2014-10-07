#!/usr/bin/env python2
'''
Choose particles positions for tracks.

Usage:
  pmovie_select.py [options] <input> <output>
  pmovie_select.py [options] --count-only <input>
  pmovie_select.py [options] --histograms <input> <output>

Options:
  --pool-size PS -p  PS      Set the size of the processor pool [default: 8].
  --xilim=XILIM              Set initial X limits, using a tuple [default: (float('-inf'),float('inf'))].
  --yilim=YILIM              Set initial Y limits, using a tuple [default: (float('-inf'),float('inf'))].
  --zilim=ZILIM              Set initial Z limits, using a tuple [default: (float('-inf'),float('inf'))].
  --select=N                 Randomly select N initial particle positions [default: 100].
  --verbose -v               Allow verbose output.
'''
from docopt import docopt;
import lspreader as rd;
import numpy as np;
import random as rnd;
import cPickle;
import itertools as it;
import sys;

def logprint(s):
    print(s);
opts = docopt(__doc__,help=True);
outname = opts['<output>'];
name = opts['<input>'];
ps = int(opts['--pool-size']);
xlim = eval(opts['--xilim']);
ylim = eval(opts['--yilim']);
zlim = eval(opts['--zilim']);
select = int(opts['--select']);
if opts['--verbose']:
    vprint = logprint;
else:
    vprint = lambda s: None;
vprint('reading in {}'.format(name));
with rd.LspOutput(name,verbose=True) as f:
    d=f._getmovie();
frame = d[0];
del d;
d=frame['data'];
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
    plt.hist(d['xi'],bins=100);
    plt.subplot(312);
    plt.title('yi');
    plt.hist(d['yi'],bins=100);
    plt.subplot(313);
    plt.title('zi');
    plt.hist(d['zi'],bins=100);
    plt.savefig(outname);
    exit(0);
out = {};

vprint('beginning selection');
vprint('xlim:{}; ylim:{}; zlim:{};'.format(xlim,ylim,zlim));
sys.stdout.flush();

goods = (xlim[0] <= d['xi']) & (d['xi'] <= xlim[1]) & \
        (ylim[0] <= d['yi']) & (d['yi'] <= ylim[1]) & \
        (zlim[0] <= d['zi']) & (d['zi'] <= zlim[1]);

vprint('we have {} matches.'.format(goods.astype(int).sum()));
vprint('selecting {}'.format(select));
for k in ['xi','yi','zi']:
    out[k] = d[k][goods];
del frame,d;
keys = [k for k in out];
d = np.array([out[k] for k in out]).T;
np.random.shuffle(d);
d=(d[:select]).T
for k,di in it.izip(keys,d):
    out[k]=di;
out['xlim']=xlim;
out['ylim']=ylim;
out['zlim']=zlim;
#making indices
out['i'] = np.arange(len(out['xi']));
vprint('outputting to {}'.format(outname));
with open(outname,'wb') as f:
    cPickle.dump(out,f,2);
pass;

