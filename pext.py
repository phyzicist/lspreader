#!/usr/bin/env python2
'''
Make generate the pext.

Usage:
  pext.py <output-pickle> (<input>...)
'''
import lspreader as rd;
import cPickle;
import numpy as np;
from docopt import docopt;
massE = 0.511e6;
def main():
    opts = docopt(__doc__,help=True);
    outname = opts['<output-pickle>']
    names = opts['<input>'];
    d=[]
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
    for k in d:
        d[k] = np.array(d[k]);
    print('calculating kinetic energy');
    d['KE'] = (np.sqrt(d['ux']**2+d['uy']**2+d['uz']**2+1)-1)*massE;
    r = np.sqrt(d['x']**2+d['y']**2+d['z']**2);    
    print('calculating polar angle ("zenith")');
    d['theta'] = np.arccos(d['z']/r);
    print('calculating azimuth');
    d['phi'] = np.arccos(d['x']/(r*np.sin(d['theta'])));
    for i,y in enumerate(d['y']):
        if y < 0.0:
            d['phi'][i] = 2*np.pi - d['phi'][i];
        pass;
    print('outputting');  
    with open(outname,"wb") as f:
        cPickle.dump(d,f,2);
    print('scatter plot in 3d with the color function "q"');
    pass;

if __name__ == '__main__':
    main();
