#!/usr/bin/env python2
'''
Parse a lsp fields or scalar file and output it as an hd5 file.

There are currently no known limitations of this script.

Usage:
  flds.py [options] <input> <output>
  flds.py [options] <input> <output> <var>...

Options:
  --help -h              Show this help.
  --verbose -v           Turn on verbosity.
  --sort=SFILE -s SFILE  Sort using the given index file.
  --hdf -H               Output to hdf instead of pickle.
  --zip -z               Compress for hdf5.
'''

import lspreader2 as rd;
from docopt import docopt;
from time import time;
import numpy as np;
from misc import h5w, mkvprint;
opts=docopt(__doc__,help=True);

vprint=mkvprint(opts);

if len(opts['<var>']) == 0:
    opts['<var>'] = False;
b=time();
d=rd.read(opts['<input>'],var=opts['<var>']);
vprint("time to read file {}: {}".format(opts['<input>'],time()-b));
vprint("read: {}".format(",".join(d.keys())));
if opts['--sort']:
    vprint("sorting using {}".format(opts['--sort']));
    sortargs = np.load(opts['--sort']);
    for k in d:
        d[k] = d[k][sortargs];
vprint("outputting to {}".format(opts['<output>']));
if opts['--hdf']:
    h5w(opts['<output>'], d,
        compression='lzf' if opts['--zip'] else None);
else:
    dump_pickle(opts['<output>'], d);
