#!/usr/bin/env python2
'''
Extract pmovie particles.

Usage:
  pmovie.py <input> <output>
'''
import lspreader2 as rd;
import cPickle as pickle
from docopt import docopt;

def _print(s):
    print(s);
vprint = _print if opts['--verbose'] else lambda s:None;

opts = docopt(__doc__,help=True);
d=rd.read(opts['<input>']);
with open(opts['<output>'],'w') as f:
    pickle.dump(d,f,2);

