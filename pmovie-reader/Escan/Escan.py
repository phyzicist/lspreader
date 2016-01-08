#!/usr/bin/env python2
'''
Scan a pmovie file for trajectories. Use this as a template.

Usage:
    ./Escan.py [options] <input> <output>

Options:
    --help -h                 Print this help.
'''

import numpy as np;
from misc import readfile;
from functools import reduce;
from docopt import docopt;

opts=docopt(__doc__,help=True);
minE = 0.120;

d=readfile(opts["<input>"],dumpfull=True);
data=frame['data'];
E   = (np.sqrt(data['ux']**2+data['uy']**2+data['uz']**2+1)-1)*0.511;
hs  = data['hash'][np.where(E>minE)];
np.save(opts['<output>'], hs);
