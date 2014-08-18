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
'''
import lspreader as rd;
import cPickle;
import numpy as np;
import itertools as itools;
from docopt import docopt;
massE = 0.511e6;

def calculate2d(x,y,d):
    print('calculating kinetic energy');
    r = np.sqrt(d['u'+x]**2+d['u'+y]**2);    
    d['KE'] = (np.sqrt(r**2+1)-1)*massE;
    print('calculating azimuth');
    d['phi'] = np.arccos(-d['u'+x]/r);
    for i,zi in enumerate(d['u'+y]):
        if zi < 0.0:
            d['phi'][i] = -d['phi'][i];
        pass;
    return d;
def calculate3d(d):
    print('calculating kinetic energy');
    r = np.sqrt(d['ux']**2+d['uy']**2+d['uz']**2);    
    d['KE'] = (np.sqrt(r**2+1)-1)*massE;
    print('calculating polar angle ("zenith")');
    d['theta'] = np.arccos(d['uz']/r);
    print('calculating azimuth');
    d['phi'] = np.arccos(-d['ux']/(r*np.sin(d['theta'])));
    for i,z in enumerate(d['uz']):
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
    latetime = float(opts['--late-time']) if opts['--late-time'] else None;
    if num_of_coords==0:
        num_of_coords=3;
    d=[];
    for name in names:
        print('reading in {}'.format(name));
        with rd.LspOutput(name) as f:
            d.append(f.get_data());
    print('stringing together');
    data=d[0];
    for i in d[1:]:
        for k in data:
            data[k].extend(i[k]);
    d=data;
    del data;
    print('cutting out duplicate times');
    #ugly c-like loop
    i = 0;
    while i<len(d['t'])-1:
        if d['t'][i] > d['t'][i+1]:
            cuti = d['t'][i+1:].index(d['t'][i])+i+2;
            for k in d:
                del d[k][i+1:cuti];
        i+=1;
    if latetime:
        print('cutting out times greater than {}'.format(latetime));
        for i,t in enumerate(d['t']):
            if t > latetime:
                for k in d:
                    del d[k][i:];
                break;
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
