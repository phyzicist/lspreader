#!/usr/bin/env python2
'''
Generate random indices for use with p4toS.py.

Usage:
  ./rndi.py <number> <input> <output>
  ./rndi.py --count-only <input>
'''

import lspreader as rd;
import random;
import sys;
import cPickle;
from docopt import docopt;
opts=docopt(__doc__,help=True);
if opts['--count-only']:
    count=True;
else:
    outname=opts['<output>'];
    take=int(opts['<number>']);
    count=False;
inname=opts['<input>'];
print('counting points');
with rd.LspOutput(inname) as f:
    if f.header['dump_type'] == 2 or f.header['dump_type'] == 3:
        if f.header['dump_type'] == 2:
            size=3;
        else:
            size=1;
        total=0;
        qs = [i[0] for i in f.header['quantities']];
        for k in range(f.header['domains']):
            f.seek(12,1);#skip *R
            nI = f.get_int();
            f.seek(4*nI,1);
            nJ = f.get_int();
            f.seek(4*nJ,1);
            nK = f.get_int();
            f.seek(4*nK,1);
            nAll=nI*nJ*nK;
            total+=nAll;
            f.seek(nAll*size*4*len(qs),1);
    else:
        raise ValueError('incorrect file type {}'.format(f.header['type']));
    pass;
print('got {}'.format(total));
if count:
    quit();
print('generating {} random numbers'.format(take));
out = [random.randint(0,total) for i in range(take)];
out.sort();
print('outputting');
with open(outname,'wb') as f:
    cPickle.dump(out,f,2);
pass;
