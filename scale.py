#!/usr/bin/env python2
'''
Scale a scalar file.

Usage:
  p4toH.py [options] <input> <scale>

Options:
  --help -h               Show this help.
  --divide -d             Divide instead of multiply.
'''

import cPickle;
import numpy as np;
from docopt import docopt;
import random;
from misc import conv;

def main():
    opts=docopt(__doc__,help=True);
    name = opts['<input>'];
    scale = conv(opts['<scale>'],default=1.0,func=float);
    with open(name,'r') as f:
        S=cPickle.load(f);
    if opts['--divide']:
        S/=scale;
    else:
        S*=scale;
    with open(name,'wb') as f:
        cPickle.dump(S,f,2);
    pass;

if __name__ == '__main__':
    main();
