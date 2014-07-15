#!/usr/bin/env python2

import lspreader as rd;
import cPickle;
import numpy as np;
import scipy.interpolate as interpol;
from docopt import docopt;
import random;


def logprint(s):
    print("{}: {}".format(name,s));

def interpolate_scalar(x,y,z,s,res=100j):
    '''Interpolates the scalar s on the grid x,y,z
       with the resolution res.'''
    Z,Y,X = np.mgrid[ min(z) : max(z) : res,
                      min(y) : max(y) : res,
                      min(x) : max(x) : res];
    S = interpol.griddata((x,y,z),s,(X,Y,Z));
    return S;

def main():
    usage="usage: ./p4toS.py <var> <rnd-pickle> <input> <output-pickle>"
    opts=docopt(usage,help=True);
    var=opts['<var>'];
    rndname=opts['<rnd-pickle>'];
    inname=opts['<input>'];
    global name;
    name = inname;
    outname=opts['<output-pickle>'];
    logprint('reading in {}'.format(inname));
    with rd.LspOutput(inname) as f:
        d = f.get_data(var=[var],pool_size=24,lazy=False);
    logprint("preparing to prune, reading in {}".format(rndname));
    with open(rndname,'rb') as f:
        rnd = cPickle.load(f);
    logprint("pruning...");
    for k in d:
        d[k] = [d[k][i] for i in rnd];
    logprint('making arrays for interpolating scalar field {}'.format(var));
    x=np.array(d['x'])*10000;
    y=np.array(d['y'])*10000;
    z=np.array(d['z'])*10000;
    s=np.log10(np.array(d[var])+0.1)#to avoid 0.0
    logprint('interpolating');
    S = interpolate_scalar(x,y,z,s);
    del x,y,z,s,d;
    logprint("dumping");
    with open(outname,"wb") as f:
        cPickle.dump(S,f,2);
    pass;

if __name__ == '__main__':
    main();
