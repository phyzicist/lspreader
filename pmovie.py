#!/usr/bin/env python2
'''
Make pmovie particles.

Usage:
  pmovie.py [options] <select-file> <input>

Options:
  --parallel=PS -p PS       Set parallel with pool size. [default: sequential]
  --subdivisions=S -s S     Run sequentially with the following subdivisions S of the pmovie. [default: 8]
  --verbose                 Run verbosely
'''
import lspreader as rd;
import cPickle;
import numpy as np;
from docopt import docopt;
import random;
Callable = rd.Callable;
import itertools as itools;
from time import time;
from misc import conv;
import multiprocessing;

def getmatches(pair):
    b,e = pair;
    print('working pair {}'.format(pair));
    bt = time();
    data = itools.izip(
        range(b,e),
        frame['xi'][b:e],
        frame['yi'][b:e],
        frame['zi'][b:e]);
    seldata = zip(
        select['i'],
        select['xi'],
        select['yi'],
        select['zi']
    );
    keys = [k for k in frame if k not in ['t','pnum','step']];
    out = {k:[] for k in keys}
    out['i'] = [];
    for i,xi,yi,zi in data:
        for j,sxi,syi,szi in seldata:
            if eq(xi,sxi) and eq(yi,syi) and eq(zi,szi):
                for k in keys:
                    out[k].append(frame[k][i]);
                out['i'].append(j);
    print('elapsed time: {}'.format(time()-bt));
    return out;

def eq(x,y):
    return abs(x-y) <= 1e-16;

opts = docopt(__doc__,help=True);
inname = opts['<input>'];
selectfile = opts['<select-file>'];
if opts['--parallel']:
    ps = 1;
    ss = int(opts['--parallel']);
else:
    ps = conv(opts['--subdivisions'],default=8,func=int);
    ss = 1;
    #reading in select file.
print('reading in select file {}'.format(selectfile));
with open(selectfile,'r') as f:
    s = cPickle.load(f);
for k in s:
    s[k]=np.array(s[k]);
select = s;#redeclaration in gloabal scope.
print('reading in {}'.format(inname));
with rd.LspOutput(inname,verbose=opts['--verbose']) as f:
    d=f._getmovie(ps);
frame = d[0];
del d[:1];
N = frame['pnum'];
l = range(0,N,N/ss)[:-1]+[N] if ss != 1 else (0,N);
ll = zip(l,l[1:]);#try this if you are uncertain.
if opts['--parallel']:
    #finally, selecting
    pool = multiprocessing.Pool(ss);
    results = pool.map(getmatches,ll);
    pool.close();
    print('stringing together');
    out = results[0];
    for i in results[1:]:
        for k in out:
            out[k].extend(i[k]);
    del results;
else:
    results = map(getmatches,ll);
    out = results[0];
#creating output name
outname='points{}.pt'.format(frame['step']);
print('outputting {}'.format(outname));
with open(outname,'w') as f:
    cPickle.dump(out,f,2);
pass;




