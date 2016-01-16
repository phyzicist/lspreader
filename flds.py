#!/usr/bin/env python2
'''
Parse a lsp fields or scalar file.

Usage:
  flds.py [options] <input> <output>
  flds.py [options] <input> <output> <var>...

Options:
  --help -h              Show this help.
  --verbose -v           Turn on verbosity.
  --npz -n               Use compressed npz over pickle.
  --reshape -r           Use rectangular reshaping for contiguous simulations.
'''
from time import time;
import numpy as np;

def rect_flds(d):
    dims = ['xs', 'ys', 'zs'];
    labels = [ key for key in d.keys()
               if key not in dims ];
    shape = [ len( np.unique(d[l]) )
              for l in dims ];
    for l in labels:
        d[l].reshape(shape);
    for l in dims:
        del d[l];
    return d;

if __name__ == "__main__":
    from docopt import docopt;
    from misc import dump_pickle;
    import lspreader as rd;
    
    vprint = mkvprint(opts);
    opts=docopt(__doc__,help=True);
    if len(opts['<var>']) == 0:
        opts['<var>'] = False;
    b=time();
    d=rd.read(opts['<input>'],var=opts['<var>'],vprint=vprint);
    vprint("time to read file {}: {}".format(opts['<input>'],time()-b));
    vprint("read: {}".format(",".join(d.keys())));
    if opts['--reshape']:
        d = rect_flds(d);
    if opts['--npz']:
        dims = ['xs', 'ys', 'zs'];
        labels = [ key for key in d.keys()
                 if key not in dims ];
        dtype = [ (l,'f4') for l in labels ];
        dim = np.array(
            [d[l] for l in dims],
            dtype=[ (l, 'f4') for l in dims ]);
        d = np.array(
            [d[l] for l in labels],
            dtype=dtype);
        np.savez_compressed(opts['<output>'],
                            dims=dim,
                            data=d);
    else:
        dump_pickle(opts['<output>'],d);
