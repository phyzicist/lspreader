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
  --X -x                  Use X.
  --Y -y                  Use Y.
  --Z -z                  Use Z.
'''
import lspreader as rd;
import cPickle;
import numpy as np;
import itertools as itools;
from docopt import docopt;
massE = 0.511e6;

def calculate2d(x,y,d):
    print('calculating kinetic energy');
    d['KE'] = (np.sqrt(d['u'+x]**2+d['u'+y]**2+1)-1)*massE;
    r = np.sqrt(d[x]**2+d[y]**2);    
    print('calculating azimuth');
    d['phi'] = np.arccos(-d[x]/r);
    for i,zi in enumerate(d[y]):
        if zi < 0.0:
            d['phi'][i] = -d['phi'][i];
        pass;
    return d;
def calculate3d(d):
    print('calculating kinetic energy');
    d['KE'] = (np.sqrt(d['ux']**2+d['uy']**2+d['uz']**2+1)-1)*massE;
    r = np.sqrt(d['x']**2+d['y']**2+d['z']**2);    
    print('calculating polar angle ("zenith")');
    d['theta'] = np.arccos(d['z']/r);
    print('calculating azimuth');
    d['phi'] = np.arccos(-d['x']/(r*np.sin(d['theta'])));
    for i,z in enumerate(d['z']):
        if z < 0.0:
            d['phi'][i] = -d['phi'][i];
        pass;
    return d;
    
def main():
    opts = docopt(__doc__,help=True);
    outname = opts['<output>']
    names = opts['<names>'];
    coords = {'x':opts['--X'],'y':opts['--Y'],'z':opts['--Z']};
    num_of_coords = len([k for k in coords if coords[k]]);
    if num_of_coords==0:
        num_of_coords=3;
    d=[];
    for name in names:
        print('reading in {}'.format(name));
        with rd.LspOutput(name) as f:
            d.append(f.get_data());
    print('cutting out duplicate times');
    for first,second in itools.izip(d[::2],d[1::2]):
        try:
            tp = first['t'][-1];
        except IndexError:
            continue;
        pass;
        #this assumes that the times are saved sequentially
        try:
            i=second['t'].index(tp);
        except ValueError:
            i=0;
        #cutting off overrun
        for k in second:
            second[k]=second[k][i:];
    data=d[0];
    print('stringing together and removing');
    for i in d[1:]:
        for k in data:
            data[k].extend(i[k]);
    d=data;
    for k in d:
        d[k] = np.array(d[k]);
    if num_of_coords == 2:
        #notice this only does x-y, y-z, and inverts the first.
        x = 'x' if opts['--X'] else 'y';
        y = 'y' if opts['--Y'] else 'z';
        d = calculate2d(x,y,d);
    elif num_of_coords ==3:
        d = calculate3d(d);
    else:
        raise RuntimeError, "derp";
    print('outputting');  
    with open(outname,"wb") as f:
        cPickle.dump(d,f,2);
    print('scatter plot in 3d with the color function "q"');
    pass;

if __name__ == '__main__':
    main();
