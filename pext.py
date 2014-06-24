#!/usr/bin/env python2

import lspreader as rd;
import cPickle;
import numpy as np;
import sys;

def main():
    usage="usage: ./p4toS.py <output-pickle> <input>..."
    if len(sys.argv) <= 2:
        print(usage);
        exit();
    outname=sys.argv[1];
    names = sys.argv[2:];
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
    d['KE'] = d['ux']**2+d['uy']**2+d['uz']**2;
    
    print('calculating azimuth');
    d['phi'] = np.arctan(d['y'] / d['x']);
    
    print('calculating polar angle ("zenith")');
    r = np.sqrt(d['x']**2+d['y']**2+d['z']**2);
    d['theta'] = d['z']/r;
    print('calculating scatter plot');
    #this implicitly deletes the positions.
    d['x'] = np.array(d['KE'] * np.sin(d['theta']) * np.cos(d['phi']));
    d['y'] = np.array(d['KE'] * np.sin(d['theta']) * np.sin(d['phi']));
    d['z'] = np.array(d['KE'] * np.cos(d['theta']));
    del d['theta'], d['phi'];
    print('outputting');  
    with open(outname,"wb") as f:
        cPickle.dump(d,f,2);
    print('scatter plot in 3d with the color function "q"');
    pass;

if __name__ == '__main__':
    main();
