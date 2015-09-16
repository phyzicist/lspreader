#!/usr/bin/env python2
'''
Parse a lsp fields or scalar file and output it as an hd5 file.

There are currently no known limitations of this script.

Usage:
  sclr.py [options] <input> <output>
  sclr.py [options] <input> <output> <var>...

Options:
  -h --help               Show this help.
  -v --verbose            Turn on verbosity.
  -x --X                  Use X in interpolation.
  -y --Y                  Use Y in interpolation.
  -z --Z                  Use Z in interpolation.
  --sort -S               Use the experimental sorting algorithm.
'''

import lspreader2 as rd;
from docopt import docopt;
import h5py as h5;
from time import time;
opts=docopt(__doc__,help=True);
def _print(s):
    print(s);
vprint = _print if opts['--verbose'] else lambda s:None;
if len(opts['<var>']) == 0:
    opts['<var>'] = False;
b=time();
d=rd.read(opts['<input>'],var=opts['<var>']);
vprint("time to read file {}: {}".format(time()-b));

if opts['--sort']:
    vprint("sorted per the passed option");
    b = time();
    sorted = np.lexsort((d['x'],d['y'],d['z']));
    for k in d:
        d[k] = d[k][sorted];
    vprint("time to sort: {}".format(time()-b));
vprint("outputting to {}".format(opts['<output>']));

with h5.File(opts['<output>'],'a') as f:
    for k in d:
        if k not in f.keys():
            f.create_dataset(k);
        f[k] = d[k];

