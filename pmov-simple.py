#!/usr/bin/env python2
'''
Extract pmovie particles.

Usage:
  pmovie.py <input> <output>
'''
import lspreader as rd;
import cPickle as pickle
from docopt import docopt;

opts = docopt(__doc__,help=True);
with rd.LspOutput(opts['<input>']) as f:
    d=f._getmovie();
with open(opts['<output>'],'w') as f:
    pickle.dump(d,f,2);

