#!/usr/bin/env python2
'''
Parse a lsp fields or scalar file and output a 3D array sampling the
given variables. The scheme used here is reading the data in, either
interpolating the scalar value (using the "nearest" point in actual data)
to a have a lower number of samples of the actual data, or performing a
N-D histogram of the data using the value of the scalar field as a histogram
weight, and then dividing by the number of points to obtain an average
over each point. The default is "nearest" interpolation.
The input to this script should be a .p4 output from lsp.

There are currently no known limitations of this script. It is sequential.

Usage:
  sclr.py [options] <input> (<var> <output>)...
  sclr.py --list-var <input>

Options:
  -h --help               Show this help.
  -v --verbose            Turn on verbosity.
  -x --X                  Use X in interpolation.
  -y --Y                  Use Y in interpolation.
  -z --Z                  Use Z in interpolation.
  --list-var              Just list the variable names.
  --xres=XRES             Set the resolution along the x direction [default: 100].
  --yres=XRES             Set the resolution along the y direction [default: 100].
  --zres=XRES             Set the resolution along the z direction [default: 100].
  --use-sort              Use the experimental sorting algorithm.
  --histogram -H          Histogram instead.
  --permute -p            Swap the order of axes for 2D data.
'''

import lspreader_old as rd;
import cPickle;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.interpolate as interp;
from docopt import docopt;
from time import time;
print("!!!WARNING: This script will be deprecated soon, within the next week.");

def logprint(s):
    global verbose;
    if verbose:
        print("{}: {}".format(name,s));
    pass;

def histogram_scalar_3d(x,y,z,s,
                        xres=100,yres=100,zres=100):
    '''Histograms the scalar s on the grid x,y,z
       with the resolutions *res along each axis.'''
    maxs = dict(x=x.max(),y=y.max(),z=z.max());
    mins = dict(x=x.min(),y=y.min(),z=z.min());
    xbins = np.linspace(mins['x'],maxs['x'],xres+1);
    ybins = np.linspace(mins['y'],maxs['y'],yres+1);
    zbins = np.linspace(mins['z'],maxs['z'],zres+1);
    t=time();
    H,_= np.histogramdd([x,y,z],bins=[xbins,ybins,zbins],weights=s);
    logprint('It took {} seconds'.format(time()-t));
    #divide to get the average
    norm=float(len(s))/float(xres*yres*zres);
    logprint('we have a norm of {}'.format(norm));
    #making bounds
    xb,yb,zb = map(lambda s: (mins[s],maxs[s]), ['x','y','z']);
    return H/norm, xb, yb, zb;

def interpolate_scalar_3d(x,y,z,s,
                          xres=100,yres=100,zres=100):
    '''Interpolates the scalar s on the grid x,y
       with the resolutions *res along each axis.
       Uses "nearest" interpolation'''
    maxs = dict(x=x.max(),y=y.max(),z=z.max());
    mins = dict(x=x.min(),y=y.min(),z=z.min());
    
    X,Y,Z = np.mgrid[ mins['x']:maxs['x']:xres*1j,
                      mins['y']:maxs['y']:yres*1j,
                      mins['z']:maxs['z']:zres*1j ];
    t=time();
    H = interp.griddata((x,y,z),s,(X,Y,Z),method='nearest');
    logprint('It took {} seconds'.format(time()-t));
    xb,yb,zb = map(lambda s: (mins[s],maxs[s]), ['x','y','z']);
    return H, xb, yb, zb;


def histogram_scalar_2d(x,y,s,
                        xres=100,yres=100):
    '''Histograms the scalar s on the grid x,y
       with the resolutions *res along each axis.'''
    maxs = dict(x=x.max(),y=y.max());
    mins = dict(x=x.min(),y=y.min());

    xbins = np.linspace(mins['x'],maxs['x'],xres+1);
    ybins = np.linspace(mins['y'],maxs['y'],yres+1);
    H,_,_ = np.histogram2d(x,y,bins=(xbins,ybins),weights=s);
    norm=float(len(s))/float(xres*yres);
    logprint('we have a norm of {}'.format(norm));
    xb,yb = map(lambda s: (mins[s],maxs[s]), ['x','y']);
    return H/norm, xb, yb;

def interpolate_scalar_2d(x,y,s,
                          xres=100,yres=100):
    '''Interpolates the scalar s on the grid x,y
       with the resolutions *res along each axis.
       Uses "nearest" interpolation'''
    maxs = dict(x=x.max(),y=y.max());
    mins = dict(x=x.min(),y=y.min());
    #ordering is to keep matlab semantics
    Y,X = np.mgrid[ mins['y']:maxs['y']:yres*1j,
                    mins['x']:maxs['x']:xres*1j];
    H = interp.griddata((y,x),s,(Y,X),method='nearest');
    xb,yb = map(lambda s: (mins[s],maxs[s]), ['x','y']);
    return H, xb, yb

def interpolate_scalar_1d(x,s,
                          res=100):
    '''Interpolates the scalar s on the grid x
       with the resolutions *res along each axis.
       Uses "nearest" interpolation'''
    maxx = x.max(); minx = x.min()
    X = np.mgrid[ minx:maxx:res*1j];
    H = interp.griddata(x,s,X,method='nearest');
    return H, (minx,maxx);

def histogram_scalar_1d(x,s,res=100):
    '''Histograms the scalar s on the grid x
       with the resolution res.'''
    xbins = np.linspace(min(x),max(x),res+1);
    H,_ = np.histogram(x,bins=xbins,weights=s);
    norm=float(len(s))/float(res);
    logprint('we have a norm of {}'.format(norm));
    return H/norm;

def main():
    global name,verbose;
    #procesing arguments
    opts=docopt(__doc__,help=True);
    var=opts['<var>'];
    outnames= opts['<output>'];
    vopairs = zip(var,outnames);
    name = opts['<input>'];

    verbose = opts['--verbose'];
    if opts['--list-var']:
        with rd.LspOutput(name, verbose=verbose, prefix=name) as f:
            names = zip(*f.header['quantities'])[0];
        for name in names:
            print(name);
        return;
    use = [];
    res = [];
    if opts['--X']:
        use.append('x');
        res.append(int(opts['--xres']));
    if opts['--Y']:
        use.append('y');
        res.append(int(opts['--yres']));
    if opts['--Z']:
        use.append('z');
        res.append(int(opts['--zres']));
    if use == []:
        use = ['x','y','z'];
        res = map(lambda k: int(opts[k]),['--xres','--yres','--zres']);
    # A couple of things to note; written in this way, whatever
    # this list (and thus, what is read) becomes, it is ordered
    # alphabetically. This is important, as this determines what
    # each resulting row and column and breadth in the output
    # array corresponds to from the actual simulation.
    #
    # It is probably worth mentioning that the xz in simulation
    # axes will be [0,1] in numpy axes, that is, it will be left-handed.
    # Using xz leads to this anyway, but it's worth reminding the reader.
    # To permute in 2D, use the --permute flag.
    logprint('reading in {}'.format(name));
    with rd.LspOutput(name, verbose=verbose, prefix=name) as f:
        d = f.get_data(var=var);
    # removing stuff that won't be read;
    vs = [v[0] for v in vopairs] +  ['x','y','z'];
    d = {k:d[k] for k in d if k in vs};

    if opts['--use-sort']:
        # making numpy array for sorting
        t = time();
        l = len(d[d.keys()[0]]);
        dt = zip(d.keys(),['f32']*len(d));
        #magic
        d = np.array(zip(*[iter(np.array(d.values()).T.ravel())]*len(d)),dtype=dt);
        d.sort(order=['x','y','z'])
        logprint("sorting took {} seconds.".format(time()-t));
        
    for v,outname in vopairs:
        logprint('histogram' if opts['--histogram'] else 'interpolating');
        o={};
        if len(use) == 3:
                f  = histogram_scalar_3d if opts['--histogram'] else interpolate_scalar_3d
                o['s'], o['x'],o['y'],o['z'] = f(d['x'],d['y'],d['z'],d[v],
                                                xres=res[0],yres=res[1],zres=res[2]);
        elif len(use) == 2:
            if opts['--permute']:
                t = use[0]; use[0] = use[1]; use[1] = t;
            f  = histogram_scalar_2d if opts['--histogram'] else interpolate_scalar_2d
            o['s'], o[use[0]], o[use[1]] = f(d[use[0]],d[use[1]],d[v],xres=res[0],yres=res[1]);
            #to remember the order
            o['0th'] = use[0]; o['1st'] = use[1];
        else:
            f  = histogram_scalar_1d if opts['--histogram'] else interpolate_scalar_1d
            o['s'], o[use[0]] = f(d[use[0]],d[v],res=res[0]);
        logprint("outputting");
        with open(outname,"wb") as f:
            cPickle.dump(o,f,2);
    pass;

if __name__ == '__main__':
    main();
