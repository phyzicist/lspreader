#!/usr/bin/env python2
'''
Sort pmovie files in place. Currently only works on the first record.

Usage:
  pmov-sort.py [options] <input>

Options:
  --help   -h            You know.
  --X=XLIM -x XLIM       Use X with limits. [default: (-30e-4, 5e-4)]
  --Y=YLIM -y YLIM       Use Y with limits. [default: None]
  --Z=ZLIM -z ZLIM       Use Z with limits. [default: (-20e-4,20e-4)]
'''
import lspreader as rd;
import cPickle as pickle
from docopt import docopt;
import numpy as np;
import misc as m;
opts=docopt(__doc__,help=True);
lims = [eval(opts[s]) for s in ["--X","--Y","--Z"]];
if lims == [None,None,None]:
    quit();
dims = [name for i,name in zip(lims,['xi','yi','zi']) if i is not None];
lims  = [lim for lim in lims if lim is not None];
f = m.readfile(opts['<input>'], dumpfull=True)
d = f[0]['data'];
keys=np.array([d[s] for s in dims]);
for i,lim in enumerate(lims):
    keys[i] -= lim[0];
sortedargs = np.lexsort(keys);
f[0]['data'] = d[sortedargs];
with open(opts['<input>'],'w') as outf:
    pickle.dump(f,outf,2);
pass
