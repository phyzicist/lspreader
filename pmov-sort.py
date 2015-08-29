#!/usr/bin/env python2
'''
Sort pmovie files in place. Currently only works on the first record. By default,
works on x and z.

Usage:
  pmov-sort.py [options] <input>

Options:
  --help -h       You know.
  --X -x          Use X. 
  --Y -y          Use Y.
  --Z -z          Use Z.
'''
import lspreader as rd;
import cPickle as pickle
from docopt import docopt;
import numpy as np;
import misc as m;
from itertools import izip;
opts=docopt(__doc__,help=True);
got = [opts[s] for s in ["--X","--Y","--Z"]];
dims=[label for s,label in zip(got,['xi','yi','zi']) if s];
if dims == []:
    dims = ['xi','zi'];
f = m.readfile(opts['<input>'], dumpfull=True)
d = f[0]['data'];
#first, prune and take only unique particles
_,uniq = np.unique(d[dims],return_index=True);
d=d[uniq];
keys=np.array([d[s] for s in dims]);
sortedargs = np.lexsort(keys);
f[0]['data'] = d[sortedargs];
#with open(opts['<input>'],'w') as outf:
#    pickle.dump(f,outf,2);
