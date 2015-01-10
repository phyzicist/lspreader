'''
Miscellaneous definitions.
'''

import cPickle as pickle;
import numpy as np;
from matplotlib import colors;

def conv(arg,default=None,func=None):
    if func:
        return func(arg) if arg else default;
    else:
        return arg if arg else default;

def read(filename,dictlabel='s', dumpfull=False):
    with open(filename,'r') as f:
        d=pickle.load(f);
    if type(d) == np.ndarray or dumpfull:
        return d;
    elif type(d) == dict:
        return d[dictlabel];
    else:
        s = str(type(d));
        errstr='Unknown pickle type "{}" loaded from file "{}".'.format(s,filename);
        raise IOError(errstr);

def cmap():
    hsv = np.array([[[0.00, 0.6, 1.0],
                     [0.10, 0.6, 1.0],
                     [0.20, 0.6, 1.0],
                     [0.30, 0.6, 1.0],
                     [0.40, 0.6, 1.0],
                     [0.50, 0.6, 1.0],
                     [0.60, 0.6, 1.0],
                     [0.70, 0.6, 1.0],
                     [0.80, 0.6, 1.0]]]);
    rgb = colors.hsv_to_rgb(hsv);
    def mk_component(cmp):
        l = len(cmp);
        inter = np.linspace(0.001,1.0,l);
        ret = [[i,j,j]  for i,j in zip(inter,cmp)];
        ret[0][1]=1.0;
        ret = [[0.0,1.0,1.0]]+ret;
        return tuple(ret);
    r = np.array(mk_component(rgb.T[0]));
    g = np.array(mk_component(rgb.T[1]));
    b = np.array(mk_component(rgb.T[2]));
    cd={'red':r,'green':g,'blue':b};
    return colors.LinearSegmentedColormap('cmap',cd, 1024);

pastel_rainbow = cmap();
