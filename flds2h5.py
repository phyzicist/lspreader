#!/usr/bin/env python2
'''
Parse a lsp fields or scalar file and output it as an hd5 file.

There are currently no known limitations of this script.

Usage:
  flds2h5.py [options] <input> <output>
  flds2h5.py [options] <input> <output> <var>...

Options:
  --help -h              Show this help.
  --verbose -v           Turn on verbosity.
  --lzf -l               Use lzf zipping on the files.
  --sort=SFILE -s SFILE  Sort using the given index file.
'''

import lspreader2 as rd;
from docopt import docopt;
import h5py as h5;
from time import time;
import numpy as np;

opts=docopt(__doc__,help=True);

def _print(s):
    print(s);
vprint = _print if opts['--verbose'] else lambda s:None;

if len(opts['<var>']) == 0:
    opts['<var>'] = False;
b=time();
d=rd.read(opts['<input>'],var=opts['<var>']);
vprint("time to read file {}: {}".format(opts['<input>'],time()-b));
if opts['--sort']:
    vprint("sorting using {}".format(opts['--sort']));
    sortargs = np.load(opts['--sort']);
    for k in d:
        d[k] = d[k][sortargs];

vprint("outputting to {}".format(opts['<output>']));
with h5.File(opts['<output>'],'a') as f:
    for k in d:
        if k not in f.keys():
            if opts['--lzf']:
                f.create_dataset(k,data=d[k],compression="lzf");
            else:
                f.create_dataset(k,data=d[k]);
        else:
            f[k] = d[k];
