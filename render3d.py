#!/usr/bin/env python
import numpy as np;
import cPickle;
from docopt import docopt;
import re;
import math;
usage='''Render volumetric render of a scalar field

Usage:
  render3d.py [options] IN_FORMAT OUT_FORMAT LABEL_FORMAT <lownum> <highnum>
  render3d.py [options] (--plot-single | -s) INFILE

Options:
   --min=MIN -n MIN             Set vmin in the volumetric plot to MIN.
   --max=MAX -x MAX             Set vmax in the volumetric plot to MAX.
   --plot-single -s             Just show, don't output, single scalar.
   --cutplane=Z                 Make cutplane along z at Z.
   --clabel=CLABEL              Set the colorbar label.
'''

def logprint(s):
    print(s);

def read_file(filename):
    print('loading file {}'.format(filename));
    with open(filename,'rb') as f:
        d=cPickle.load(f);
    return d;

def tozero(v):
    if math.isnan(v):
        return 0;
    else:
        return v;

def zero_nan(S):
    return np.array([[[tozero(k) for k in j] for j in i] for i in S]);

def initial_plot(mlab,fname,vlim,
                 clabel=None,cutplane=None,label=None):
    '''
    Creates the first plot.
    '''
    ret = {};
    S=read_file(fname);
    S=np.log10(zero_nan(S)+0.1);
    fig=mlab.figure(size=(1280,1024));
    ret['fig']=fig;
    fig.scene.disable_render=True;
    src=mlab.pipeline.scalar_field(S);
    #volume rendering
    v=mlab.pipeline.volume(src,vmin=vlim[0],vmax=vlim[1]);
    ret['v']=v;
    #cut plane 1
    if cutplane is not None:
        zp=mlab.pipeline.image_plane_widget(src,plane_orientation="z_axes",
                                            slice_index=cutplane,
                                            vmin=vlim[0],vmax=vlim[1]);
        ret['zp'] = zp;
        
    mlab.view(azimuth=150,elevation=45,focalpoint='auto',distance='auto');
    mlab.scalarbar(object=v,title=clabel);
    #Unfortunately, the volume renderer doesn't work out of the box.
    #The main issue is the colorbar. So, we do this instead.
    e=mlab.get_engine();
    module_manager = e.scenes[0].children[0].children[0];
    module_manager.scalar_lut_manager.use_default_range = False;
    module_manager.scalar_lut_manager.data_range = np.array([vmin, vmax]);
    if label is not None:
        t=mlab.text(0.075,0.875,label,width=0.1);
        ret['t'] = t;
    #begin rendering
    fig.scene.disable_render=False;
    return ret;

def plot(names,vlim,
         clabel=None,cutplane=None):
    print("loading mlab");
    import mayavi.mlab as mlab;
    mlab.options.offscreen = True;
    outname,inname,label = names[0];
    names = names[1:];
    print("processing first file {}...".format(inname));
    d=initial_plot(mlab,inname,vlim,clabel=clabel,
                   cutplane=cutplane,label=label);
    print("saving {}".format(outname));
    mlab.savefig(outname,size=(1280,1024));
    for outname,inname,label in names:
        print("processing {}...".format(inname))
        S=read_file(inname);
        S=np.log10(zero_nan(S)+0.1);
        d['fig'].scene.disable_render=True;
        d['v'].mlab_source.set(scalars=S);
        if cutplane:
            d['zp'].mlab_source.set(scalars=S);
        d['t'].text=label;
        d['fig'].scene.disable_render=False;
        print("saving {}".format(outname));
        mlab.savefig(outname,size=(1280,1024));
    pass



def plot_single(inname,vlim,
                clabel=None,cutplane=None):
    print("loading mlab");
    import mayavi.mlab as mlab;
    initial_plot(mlab,inname,vlim,cutplane=cutplane);
    print('plotting');
    mlab.show();

if __name__=="__main__":
    #reading in arguments
    opts=docopt(usage,help=True);
    vmin = float(opts['--min']) if opts['--min'] else -1;
    vmax = float(opts['--max']) if opts['--max'] else 23.5;
    clabel = opts['--clabel'] if opts['--clabel'] else 'log10 of electron density';
    if not opts['--plot-single']:
        if '*' not in opts['IN_FORMAT'] or '*' not in opts['OUT_FORMAT']:
            print(usage);
            quit(1);
        lownum=int(opts['<lownum>']);
        highnum=int(opts['<highnum>']);
        fmt='{{:0>{}}}'.format(len(opts['<highnum>']))
        out_fmt=re.sub(r'\*',fmt,opts['OUT_FORMAT']);
        in_fmt =re.sub(r'\*',fmt,opts['IN_FORMAT']);
        label_fmt=re.sub(r'\*',fmt,opts['LABEL_FORMAT']);
        numrange = range(lownum,highnum+1);
        outnames = [out_fmt.format(i) for i in numrange];
        innames = [in_fmt.format(i) for i in numrange];
        labels = [label_fmt.format(i) for i in numrange];
        files = zip(outnames,innames,labels);
        plot(files,(vmin,vmax),cutplane=opts['--cutplane'],clabel=clabel);
    else:
        plot_single(opts['INFILE'],(vmin,vmax),cutplane=opts['--cutplane'],clabel=clabel);
    pass;
pass;
