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

test = lambda d,k: k in d and d[k];

def readfile(filename,dictlabel='s', dumpfull=False):
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
    pass;

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

def cmap(r,g=None,b=None,a=None):
    if g is None and b is None:
        r,g,b = r;
    cd = {'red':r,'green':g,'blue':b};
    if a is not None:
        cd.update({'alpha':a});
    return colors.LinearSegmentedColormap('cmap', cd, 1024);
def list2cmap(l):
    l=np.array(l);
    
_p_h = np.linspace(0.725, 0.0, 9);
_p_s = np.ones(_p_h.shape)*0.6;
_p_v = np.ones(_p_h.shape);
_p_hsv = np.array([_p_h,_p_s,_p_v]);
_p_rgb= rgbfromhsv(_p_hsv);
_p_rgb,_p_a =whiteoutzero(_p_rgb,alpha=True);
#_p_rgb_nozero = rgbfromhsv(_p_hsv,split=False,whitezero=False);

pastel        = cmap(_p_rgb);
pastel_clear  = cmap(_p_rgb,a=_p_a);
#pastel_nozero = cmap(_pastel_rgb_nozero);
#pastel_b2r    = cmap(rgbfromhsv(
#    np.array([np.linspace(0.0,0.725,9),_pastel_s,_pastel_v])
#));
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
