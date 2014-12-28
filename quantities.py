#!/usr/bin/env python2
'''
Output the efficiency of a sim from a pext.py output.
All units are SI.

Usage:
  quantities.py [options] <input>

Options:
  --2D             Calculate 2D quantities instead.
  --E-cut=ECUT     Cutoff at this energy in MeV. [default: 0.0]
  --E_0=E0         Enter the peak electric field. [default: 4.763e12]
  --W=SPOTSIZE     Set the spotsize. [default: 2.26e-6]
  --T=PERIOD       Set the laser period. [default: 30e-15]
'''
import numpy as np;
import math as m;
import cPickle;
from docopt import docopt;
from misc import conv;

e_0 = 8.85418782e-12;
T = 30e-15;# in seconds
c = 2.99792458e8;
e = 1.60217657e-19

opts = docopt(__doc__,help=True);

E_0 = float(opts['--E_0'])
ecut = float(opts['--E-cut'])
w = float(opts['--W'])

with open(opts['<input>'],'rb') as f:
    d = cPickle.load(f);

ecut*=1e6
q,KE = zip(*[ (q,d['KE'][i]) for i,q in enumerate(d['q']) if d['KE'][i] > ecut]);
q=np.array(q);
KE=np.array(KE);
q*=1e-6;
KE  = sum(-q*KE);

print('total charge: {} pC'.format(sum(q)*1e12));
#calculated in SI...
if opts['--2D']:
    E_laser = w * m.sqrt(m.pi/2) * (c*e_0*E_0**2)/2 * T;
else:
    E_laser = w**2 * (m.pi/2) * (c*e_0*E_0**2)/2 * T;
print('efficiency is {}'.format(KE/E_laser));
print('the pulse energy is {} J'.format(E_laser));
