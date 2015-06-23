#!/usr/bin/env python2
'''
Output the efficiency of a sim from a pext.py output.
All units are SI.

Usage:
  quantities.py [options] <input>

Options:
  --2D -2                    Calculate 2D quantities instead.
  --E-cut=ECUT -e ECUT       Cutoff at this energy in MeV. [default: 0.0]
  --E_0=E0 -E E0             Enter the peak electric field in V/m. [default: 4.763e12]
  --W=SPOTSIZE -w SPOTSIZE   Set the spotsize meters. [default: 2.26e-6]
  --T=FWHM -t FWHM           Set the Full-Width Half-Max in seconds. [default: 30e-15]
  --angle=ANGLE -a ANGLE     Restrict the angle. In 2D, this is just phi; in 3D, it is solid angle.
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

E_0 = float(opts['--E_0']);

ecut = float(opts['--E-cut'])*1e6
w = float(opts['--W'])

with open(opts['<input>'],'rb') as f:
    d = cPickle.load(f);
#q's and KE's
#restricting based on energies
good = d['KE'] > ecut;
if opts['--angle']:
    angle = float(opts['--angle']);
    if opts['--2D']:
        good &= np.abs(d['phi']) > np.pi - angle/180*np.pi;
    else:
        good &= np.cos(angle/180*np.pi) < -np.sin(d['theta'])*np.cos(d['phi']);
#restricting
q  = d['q'][good];
KE = d['KE'][good];

q*=1e-6;
KE  = sum(-q*KE);

print('total charge: {} {}'.format(sum(q)*1e12,'pC/cm' if opts['--2D'] else 'pC'));
#calculated in SI...
if opts['--2D']:
    E_laser = w * m.sqrt(m.pi/2) * (c*e_0*E_0**2)/2 * T*1e-2;
else:
    E_laser = w**2 * (m.pi/2) * (c*e_0*E_0**2)/2 * T;
print('efficiency is {}'.format(KE/E_laser));
print('the pulse energy is {} J'.format(E_laser));
