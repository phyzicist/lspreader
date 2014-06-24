#!/usr/bin/env python2

import lspreader as rd;
import random;
import sys;
import cPickle;

usage="usage: ./rndi.py <number> <input> <output>"
if len(sys.argv) != 4:
    print(usage);
    exit();
take=int(sys.argv[1]);
inname=sys.argv[2];
outname=sys.argv[3];
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
print('generating {} random numbers'.format(take));
out = [random.randint(0,total) for i in range(take)];
out.sort();
print('outputting');
with open(outname,'wb') as f:
    cPickle.dump(out,f,2);
pass;
