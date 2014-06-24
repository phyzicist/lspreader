#!/usr/bin/env python2

import lspreader as rd;
import cPickle;
import numpy as np;
import scipy.interpolate as interpol;
import sys;
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
    usage="usage: ./pickletoS.py <var> <take> <input> <output-pickle>"
    if len(sys.argv) != 5:
        print(usage);
        exit();
    var=sys.argv[1];
    if sys.argv[2] == 'all':
        take=None;
    else:
        take=int(sys.argv[2]);
    inname=sys.argv[3];
    global name;
    name = inname;
    outname=sys.argv[4];
    logprint('reading in {}'.format(inname));
    with open(inname,'rb') as f:
        d = cPickle.load(f);
    if take:
        logprint("preparing to prune, zipping");
        keys = [i for i in d];
        data = [data for k,data in d.items()];
        del d;
        data = zip(*data);
        logprint('will take {} elements out of {}'.format(take,len(data)));
        logprint('sampling');
        rm=random.sample(data, take);
        logprint('unzipping');
        del data;
        dl=zip(*rm);
        d=dict(zip(keys,dl));
        del dl;
    logprint('making arrays for interpolating scalar field {}'.format(var));
    x=np.array(d['x']);
    y=np.array(d['y']);
    z=np.array(d['z']);
    s=np.array(d[var]);#to avoid 0.0
    logprint('interpolating');
    S = interpolate_scalar(x,y,z,s);
    del x,y,z,s,d;
    logprint("dumping");
    with open(outname,"wb") as f:
        cPickle.dump(S,f,2);
    pass;

if __name__ == '__main__':
    main();
