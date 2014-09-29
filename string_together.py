#!/usr/bin/env python2
'''
String points into paths.

Usage:
  string.py [options] <input-fmt> <output-fmt> <lownum> <highnum> [<step>]

Options:
  --xdim=XFACTOR         Multiply by XFACTOR. [default: (-30e-4,5e-4)]
  --ydim=YFACTOR         Multiply by YFACTOR. [default: (-20e-4,20e-4)]
  --zdim=ZFACTOR         Multiply by ZFACTOR. [default: (-20e-4,20e-4)]
  --xsize=XSIZE          The size of the scalar x axis. [default: 100]
  --ysize=YSIZE          The size of the scalar y axis. [default: 100]
  --zsize=ZSIZE          The size of the scalar z axis. [default: 100]
  --verbose -v           Be verbose.
'''

from docopt import docopt;
import cPickle;
from itertools import izip,imap;
import numpy as np;
import re;
opts=docopt(__doc__,help=True);

if '*' not in opts['<input-fmt>']:
    print(__doc__);
    quit(1);
xd=eval(opts['--xdim']);
xs=float(opts['--xsize']);
yd=eval(opts['--ydim']);
ys=float(opts['--ysize']);
zd=eval(opts['--zdim']);
zs=float(opts['--zsize']);
scale_gen = lambda d,s: lambda x: (x-d[0])*s/(d[1]-d[0]);

xconv = scale_gen(xd,xs);
yconv = scale_gen(yd,ys);
zconv = scale_gen(zd,zs);

lownum=int(opts['<lownum>']);
highnum=int(opts['<highnum>']);
numrange = range(lownum,highnum+1);
if opts['<step>']:
    step = int(opts['<step>']);
    numrange = numrange[::step];
fmt='{{:0>{}}}'.format(len(opts['<highnum>']))
in_fmt  = re.sub(r'\*',fmt,opts['<input-fmt>']);
out_fmt = re.sub(r'\*',fmt,opts['<output-fmt>']);
innames  = [in_fmt.format(i) for i in numrange];
outnames = [out_fmt.format(i) for i in numrange];
#reading first file.
print('reading first file {}'.format(innames[0]));
with open(innames[0],'r') as f:
    d = cPickle.load(f);
N = len(d['x']);

tms = highnum-lownum;
data = np.zeros((3,N,highnum-lownum));
mkpoints = lambda dd: imap(lambda e: (e[0],xconv(e[1]),yconv(e[2]),zconv(e[3])),
              izip(dd['i'],d['x'],dd['y'],dd['z']));
for i,x,y,z in mkpoints(d):
    data[0,i,0] = x;
    data[1,i,0] = y;
    data[2,i,0] = z;
with open(out_fmt.format(numrange[0]),'w') as f:
    cPickle.dump(data[:,:,:1],f,2);
for timestep,inname,outname in zip(range(tms),innames,outnames)[1:]:
    print('reading file {}'.format(inname));
    with open(inname,'r') as f:
        d = cPickle.load(f);
    for i,x,y,z in mkpoints(d):
        data[0,i,timestep] = x;
        data[1,i,timestep] = y;
        data[2,i,timestep] = z;
    print('outputting to file {}'.format(outname));
    with open(outname,'w') as f:
        cPickle.dump(data[:,:,:timestep+1],f,2);
    pass;
pass;