#!/usr/bin/env python2
import numpy as np;
import matplotlib;
matplotlib.use('Agg');
rcdict={'text.color': 'w',
        'axes.labelcolor':'w',
        'xtick.color':'w',
        'ytick.color':'w'}
matplotlib.rcParams.update(rcdict);
import yt.mods as ymods;
import matplotlib.pyplot as plt;
import cPickle;
import math;
from docopt import docopt;


def tozero(v):
    if math.isnan(v):
        return 0;
    else:
        return v;

def zero_nan(S):
    return np.array([[[tozero(k) for k in j] for j in i] for i in S]);


def read_file(filename):
    print('loading file {}'.format(filename));
    with open(filename,'rb') as f:
        d=cPickle.load(f);
    return zero_nan(d);

def mk_cam(pf,scale=15.0):
    #getting extrema
    mn,mx = pf.h.all_data().quantities['Extrema']('Density')[0];
    #make transfer function
    transfer_function = ymods.ColorTransferFunction((mn,mx),grey_opacity=True);
    #setting camera parameters
    viewing_direction = [0.005,0.005,-0.005];
    center = pf.domain_center;
    width = 1.5*pf.domain_width;
    pixels_per_side = 1024;
    north = [0,0,1];
    cam = pf.h.camera(center, viewing_direction, width, pixels_per_side,
                      transfer_function, fields=['Density'], north_vector=north,
                      steady_north=True,  no_ghost=False, sub_samples=5,
                      log_fields=False);
    cam.transfer_function.map_to_colormap(mn, mx, scale=scale,
                                          colormap='rainbow');
    return cam;
def offaxis(pf):
    viewing_direction = [0.005,0.005,-0.005];
    center = pf.domain_center;
    width = 1.5*pf.domain_width;
    pixels_per_side = 1024;
    return ymods.off_axis_projection(pf,center,viewing_direction,width,
                                     pixels_per_side, "Density",
                                     no_ghost=False);
    pass;
def output(image,mn,mx,outname,t):
    fig=plt.figure()
    plt.imshow(image,cmap='rainbow');
    c=plt.colorbar()
    c.set_label('log10 of electron density');
    c.set_ticks(np.linspace(0,1,10));
    labels = np.linspace(mn,mx,10);
    labels = ['{:.3}'.format(i) for i in labels]
    c.set_ticklabels(labels);
    cbytick_obj = plt.getp(c.ax.axes, 'yticklabels');
    plt.setp(cbytick_obj, color='w')
    plt.axis('off');
    plt.text(10,100,'t = {} fs'.format(t));
    plt.savefig(outname,facecolor='black')


usage='''Render volumetric render of a scalar field

Usage:
  ./plot3d_yt.py <infile> <outfile> <time>
'''

    
def main():
    opts=docopt(usage,help=True);
    d = dict(Density=read_file(opts['<infile>']));
    bbox = np.array([[-0.0020,0.0020],[-0.0020,0.0020],[-0.0030,0.0005]]);
    pf = ymods.load_uniform_grid(d,d['Density'].shape,1.0,bbox=bbox);
    mn,mx = pf.h.all_data().quantities['Extrema']('Density')[0];
    cam = mk_cam(pf,scale=4.0);
    image=cam.snapshot();
    output(image,mn,mx,opts['<outfile>'],opts['<time>']);
    pass
pass;
if __name__ == '__main__':
    main();
