#!/usr/bin/env python2
'''
Show conversion of energy.

Usage:
  conv.py [options] <infile>

Options:
  --help -h             Output help.
  --bins=BINS -b BINS   Number of bins. [default: 100]
  --2D                  For 2D.
'''
import math as m;
import matplotlib.pyplot as plt;
import numpy as np;
import cPickle as pickle;
import sys;
import re;
from misc import readfile
from docopt import docopt;
import itertools as it;
opts = docopt(__doc__,help=True);
d = readfile(opts['<infile>'],dumpfull=True);
#cut by azimuth
good = (d['phi'] > np.pi/2) | (d['phi'] < -np.pi/2);

s = (d['KE']*d['q']*(-1e-3))[good]; #to mililjoule
t = d['t'][good];
h,t = np.histogram(t, weights=s, bins=int(opts['--bins']));
# I'm sure there is a function for this, but my internet is out...
for i,_ in enumerate(h[1:]):
    h[i] += h[i-1];
plt.plot(t[:-1],h);
plt.ylabel("accumulated $KE$ "+ ("(mJ)" if not opts['--2D'] else "mJ/cm"));
plt.xlabel("$t$ (ns)");
plt.show();
