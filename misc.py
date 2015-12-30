'''
Miscellaneous definitions.
'''

import cPickle as pickle;
import numpy as np;
from matplotlib import colors;
import colormaps

def conv(arg,default=None,func=None):
    if func:
        return func(arg) if arg else default;
    else:
        return arg if arg else default;

test = lambda d,k: k in d and d[k];

def readfile(filename, stuff='s',
             dumpfull=False,hdf=False,
             group=None):
    if hdf:
        with h5.File(filename, "rb") as f:
            if group:
                g=f[group];
            else:
                g=f;
            if not stuff or dumpfull:
                d={l:g[l] for l in g.keys()};
            elif type(stuff) == str:
                d = g[stuff];
            else:
                d={l:g[l] for l in stuff};
    else:
        with open(filename, "rb") as f:
            d=pickle.load(f);
        if not stuff or dumpfull:
            pass
        elif type(stuff) == str:
            d = d[stuff];
        else:
            d={l:d[l] for l in stuff};
    return d;

def dump_pickle(name, obj):
    with open(name,"wb") as f:
        pickle.dump(obj,f,2);
    pass;

def chunks(l,n):
    return [l[x:x+n] for x in xrange(0, len(l), n)];

def mkvprint(opts):
    def vprint(s):
        print(s);
    if opts['--verbose']:
        return vprint;
    else:
        return lambda s: None;

#colormap stuff
cmap_min = 0.001;
def rgbfromhsv(h,s=None,v=None):
    if  s is None and v is None:
        hsv = h.T;
    elif s is not None and v is not None:
        hsv = np.array([h,s,v]).T;
    else:
        raise ValueError("Usage of rgbfromhsv is incorrect.");
    return colors.hsv_to_rgb(hsv);

def whiteoutzero(rgb,alpha=False):
    rgb=np.array([rgb,rgb]);
    rgb[0,0,:]=1.0;
    inter = np.linspace(cmap_min, 1.0, rgb.shape[1]);
    inter = np.array([inter,inter,inter]);
    rgb = np.concatenate(([inter.T],rgb));
    top = np.array([[[0.0,1.0,1.0]]]*3);
    rgb = np.concatenate((top,rgb.T),axis=1);
    if alpha:
        a = np.copy(rgb[0]);
        a[:,1:]=1.0;
        a[0,:]=0.0;
        return rgb,a;
    else:
        return rgb;
    pass;

def listtoclear(l):
    rgb,a = whiteoutzero(l,alpha=True);
    return cmap(rgb,a=a);

def cmap(r,g=None,b=None,a=None):
    if g is None and b is None:
        r,g,b = r;
    cd = {'red':r,'green':g,'blue':b};
    if a is not None:
        cd.update({'alpha':a});
    return colors.LinearSegmentedColormap('cmap', cd, 1024);
    
_p_h = np.linspace(0.725, 0.0, 9);
_p_s = np.ones(_p_h.shape)*0.6;
_p_v = np.ones(_p_h.shape);
_p_hsv = np.array([_p_h,_p_s,_p_v]);
_p_rgb= rgbfromhsv(_p_hsv);
_p_rgb,_p_a =whiteoutzero(_p_rgb,alpha=True);

_pl_rgb_r,_pl_a_r = whiteoutzero(colormaps._plasma_data[::-1],alpha=True);

_vi_rgb,_vi_a = whiteoutzero(colormaps._viridis_data,alpha=True);
_vi_rgb_r,_vi_a_r = whiteoutzero(colormaps._viridis_data[::-2],alpha=True);

pastel        = cmap(_p_rgb);
pastel_clear  = cmap(_p_rgb,a=_p_a);
plasma_clear  =  listtoclear(colormaps._plasma_data);
plasma_clear_r = listtoclear(colormaps._plasma_data[::-1]);
magma_clear  =  listtoclear(colormaps._magma_data);
magma_clear_r = listtoclear(colormaps._magma_data[::-1]);
viridis_clear =  listtoclear(colormaps._viridis_data);
viridis_clear_r =listtoclear(colormaps._viridis_data[::-1]);
def mkstrip(rgb,vmin,vmax,val,
            strip=[1.0,0.0,0.0],log10=False):
    if log10:
        val=np.log10(val);
        vmin=np.log10(vmin);
        vmax=np.log10(vmax);
    val = (val-vmin)/(vmax-vmin);
    strip = np.array(strip);
    rng = rgb[0,:,0];
    
    for i,v in enumerate(rng):
        if v > val: break;
    #finding the good point
    dx = (val-rng[i-1])/(rng[i]-rng[i-1]);
    ep = 0.03
    #interpolating
    b = (rgb[:,i,:]-rgb[:,i-1,:])*(dx-ep)+rgb[:,i-1,:];
    a = (rgb[:,i,:]-rgb[:,i-1,:])*(dx+ep)+rgb[:,i-1,:];
    #setting the strip color
    b[:,2] = a[:,1] = strip;
    c = np.concatenate([b.T,a.T]).T.reshape(3,2,3);
    #and placing it inside rgb
    return np.concatenate((rgb[:,:i,:],c,rgb[:,i:,:]),axis=1);

def mkstrip_cmap(vmin, vmax, val, strip=[1.0,0.0,0.0],log10=False):
    return cmap(
        mkstrip(_pastel_rgb,vmin, vmax, val,strip,log10=log10)
    );
