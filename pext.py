#!/usr/bin/env python2
'''
Read in pext files and output a dict. Assumes that the names
are in order of the restart files. It still isn't general enough
for any setup (which axes  you want for the angles you want
to measure phi and theta wrt). Currently, it works for angle
wrt -x for 3d and -x for 2d.

Usage:
  pext.py [options] <output> <names>...

Options:
  --X -x                    Use X.
  --Y -y                    Use Y.
  --Z -z                    Use Z.
  --late-time=TIME -l TIME  Cut out after this time.
  --reverse -r              Reverse Y and Z.
'''
import lspreader as rd;
import cPickle;
import numpy as np;
import itertools as itools;
from docopt import docopt;
massE = 0.511e6;

def calculate2d(x,y,d,reverse=False):
    r = np.sqrt(d['u'+x]**2+d['u'+y]**2);    
    d['KE'] = (np.sqrt(r**2+1)-1)*massE;
    if reverse:
        d['phi'] = np.arctan2(d['u'+y],d['u'+x]);
    else:
        d['phi'] = np.arctan2(-d['u'+x],d['u'+y]);
    return d;

def calculate3d(d,reverse=False):
    r = np.sqrt(d['ux']**2+d['uy']**2+d['uz']**2);    
    d['KE'] = (np.sqrt(r**2+1)-1)*massE;
    if reverse:
        d['theta'] = np.arccos(d['uy']/r);
        d['phi'] = np.arctan2(-d['uz'],d['ux']);
    else:
        d['theta'] = np.arccos(d['uz']/r);
        d['phi'] = np.arctan2(d['uy'],d['ux']);
    return d;
    
def main():
    opts = docopt(__doc__,help=True);
    outname = opts['<output>']
    names = opts['<names>'];
    coords = {'x':opts['--X'],'y':opts['--Y'],'z':opts['--Z']};
    num_of_coords = len([i for i in coords.values() if i]);
    latetime = float(opts['--late-time']) if opts['--late-time'] else None;
    if num_of_coords==0:
        num_of_coords=3;
    d=[];
    for name in names:
        print('reading in {}'.format(name));
        with rd.LspOutput(name) as f:
            d.append(f.get_data());
    print('cutting out duplicate times');
    ch = lambda i,j: i[ i['t'] < j['t'][0] ]
    try:
        d[:-1] = [ch(i,j) for i,j in zip(d[:-1],d[1:])]
    except IndexError:
        #this means that one of the arrays has nothing in it, ie, nothing
        #exited that plane. Ignore this.
        pass;
    #concatenating
    d = np.concatenate(d);
    if latetime:
        print('cutting out times greater than {}'.format(latetime));
        good = d['t'] <= latetime;
        d = d[good];
    #turn arrays into dict
    d = {k:d[k] for k in d.dtype.names};
    #calculating based on the number of dimensions.
    if num_of_coords == 2:
        #notice this only does x-y, y-z, and inverts the first.
        x = 'x' if opts['--X'] else 'y';
        y = 'y' if opts['--Y'] else 'z';
        d = calculate2d(x,y,d,opts['--reverse']);
    elif num_of_coords ==3:
        d = calculate3d(d,opts['--reverse']);
    else:
        raise RuntimeError, "derp";
    print('outputting');
    with open(outname,"wb") as f:
        cPickle.dump(d,f,2);
    pass;

if __name__ == '__main__':
    main();
