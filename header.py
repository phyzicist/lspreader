#!/usr/bin/env python2
'''
Just get the header of a .p4 and output it.

Usage:
  header.py [options] <input>

Options:
  --help -h               Show this help.
  --verbose v             Turn on verbosity.
  --size -s               Output the size of the header in bytes.
'''

import lspreader as rd;
import cPickle;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.interpolate as interp;
from docopt import docopt;
from time import time;


opts=docopt(__doc__,help=True);
verbose = opts['--verbose'];
name = opts['<input>'];
with rd.LspOutput(name, verbose=verbose, prefix=name) as f:
    h=f.header;
    if opts['--size']:
        print('size:{}'.format(f.file.tell()));
print(h);
