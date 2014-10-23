#!/usr/bin/env python
'''Render volumetric render of a scalar field.

Usage:
  render3d.py [options] IN_FORMAT OUT_FORMAT LABEL_FORMAT <lownum> <highnum> [<step>]
  render3d.py [options] (--plot-single | -s) INFILE

Options:
   --min=MIN -n MIN              Set vmin in the volumetric plot to MIN.
   --max=MAX -x MAX              Set vmax in the volumetric plot to MAX.
   --plot-single -s              Just show, don't output, single scalar.
   --cutplane=Z                  Make cutplane along z at Z.
   --clabel=CLABEL               Set the colorbar label. [default: log10 of electron density (cm^-3)]
   --zero-x=XRANGE               Set a range of x to zero, as a python tuple.
   --zero-y=YRANGE               Set a range of y to zero, as a python tuple.
   --zero-z=ZRANGE               Set a range of z to zero, as a python tuple.
   --azimuth=AZIMUTH             Set the azimuth [default: 135].
   --polar=POLAR                 Set the polar angle [default: 125].
   --roll=ROLL                   Roll after setting the angle [default: 0].
   --trajectories=TRAJ_FORMAT    Look for trajectories files.
'''
import numpy as np;
import cPickle;
from docopt import docopt;
import re;
import math;
from misc import conv;

def logprint(s):
    print(s);

def read_file(filename):
    print('loading file {}'.format(filename));
    with open(filename,'rb') as f:
        d=cPickle.load(f);
    return d;

def zero_range(S,zeros):
    S[zeros[0][0]:zeros[0][1],zeros[1][0]:zeros[1][1],zeros[2][0]:zeros[2][1]]=0.0;
    return S;

def initial_plot(mlab,fname,vlim,angle,**kwargs):
    '''
    Creates the first plot.
    '''
    ret = {};
    S=read_file(fname);
    S=np.nan_to_num(S);
    if 'zeros' in kwargs and kwargs['zeros'] is not None:
        S=zero_range(S,kwargs['zeros'])
    S=np.log10(S+0.1);
    fig=mlab.figure(size=(1280,1024));
    ret['fig']=fig;
    fig.scene.disable_render=True;
    src=mlab.pipeline.scalar_field(S);
    #volume rendering
    v=mlab.pipeline.volume(src,vmin=vlim[0],vmax=vlim[1]);
    #creating custom opacity transfer function
    from tvtk.util.ctf import PiecewiseFunction;
    otf = PiecewiseFunction();
    otf.add_point(vlim[0],    0.0);
    otf.add_point(vlim[0]+(vlim[1]-vlim[0])*0.9,0.06);
#    otf.add_point(vlim[0]+(vlim[1]-vlim[0])*0.99,0.1);
    otf.add_point(vlim[1],    0.5);
    v._otf = otf;
    v._volume_property.set_scalar_opacity(otf);
    ret['v']=v;
    #cut plane 1
    if 'cutplane' in kwargs:
        print('it is {}'.format(kwargs['cutplane']));
        zp=mlab.pipeline.image_plane_widget(src,plane_orientation="z_axes",
                                            slice_index=kwargs['cutplane'],
                                            vmin=vlim[0],vmax=vlim[1]);
        ret['zp'] = zp;
    mlab.view(elevation=angle[0],azimuth=angle[1],
              focalpoint='auto',distance='auto');
    
    mlab.roll(angle[2]);
    if 'clabel' in kwargs:
        mlab.scalarbar(object=v,title=kwargs['clabel']);
    else:
        mlab.scalarbar(object=v);
    #Unfortunately, the volume renderer doesn't work out of the box.
    #The main issue is the colorbar. So, we do this instead.
    e=mlab.get_engine();
    module_manager = e.scenes[0].children[0].children[0];
    module_manager.scalar_lut_manager.use_default_range = False;
    module_manager.scalar_lut_manager.data_range = np.array([vmin, vmax]);
    if 'label' in kwargs:
        t=mlab.text(0.075,0.875,kwargs['label'],width=0.1);
        ret['t'] = t;
    #show the nubs
    oa=mlab.orientation_axes();
    oa.marker.set_viewport(0,0,0.4,0.4);
    #begin rendering
    if 'render' not in kwargs or not kwargs['render']:
        fig.scene.disable_render=False;
    else:
        fig.scene.disable_render=False;
    fig.scene.camera.zoom(1.3);
    return ret;

def mk_trajectories(mlab,traj_name):
    with open(traj_name,'r') as f:
        traj = cPickle.load(f);#getting x,y,z arrays.
    #the shape of each dimension is d[x|y|z,particle,step]
    steps = len(traj[0,0,:]);
    N = len(traj[0,:,0]);
    sarray=np.arange(steps);
    conn_1 = np.array(zip(sarray[:-1],sarray[1:]));
    conn = np.vstack([conn_1+steps*i for i in range(N)])
    return (traj,conn);

def plot(names_list,vlim,angle,
         **kwargs):
    print("loading mayavi");
    import mayavi.mlab as mlab;
    mlab.options.offscreen = True;
    first = names_list[0];
    names_list = names_list[1:];
    print("processing first file {}...".format(first['in']));
    d=initial_plot(mlab,first['in'],vlim,angle,
                   label=first['label'],render='traj' in first,
                   **kwargs);
    if 'traj-files' in kwargs:
        firsttr,lasttr = kwargs['traj-files'];
        #reading the full files
        traj,con = mk_trajectories(mlab,lasttr);
        #getting starting length
        with open(firsttr,'r') as f:
            starti = len(cPickle.load(f)[0,0,:]);
        #making current part of the tracks.
        tr_cur = np.zeros(traj.shape);
        tr_cur.fill(np.nan);
        #setting up first parts
        tr_cur[:,:,:starti] = traj[:,:,:starti];
        x = np.hstack(tr_cur[0,:,:]);
        y = np.hstack(tr_cur[1,:,:]);
        z = np.hstack(tr_cur[2,:,:]);
        tracks_src=mlab.pipeline.scalar_scatter(x,y,z,
                                                np.ones(x.shape));
        tracks_src.mlab_source.dataset.lines = con;
        lines = mlab.pipeline.stripper(tracks_src);
        tracks = mlab.pipeline.surface(lines,line_width=2,opacity=1.0);
        starti+=1;
        d['fig'].scene.disable_render = False;
    print("saving {}".format(first['out']));
    mlab.savefig(first['out'],size=(1280,1024));
    for i,names in enumerate(names_list):
        print("processing {}...".format(names['in']))
        S=read_file(names['in']);
        S=np.nan_to_num(S)
        if 'zeros' in kwargs and kwargs['zeros'] is not None:
            S=zero_range(S,kwargs['zeros']);
        S=np.log10(S+0.1);
        d['fig'].scene.disable_render=True;
        d['v'].mlab_source.set(scalars=S);
        if 'cutplane' in kwargs:
            d['zp'].mlab_source.set(scalars=S);
        d['t'].text=names['label'];
        if 'traj-files' in kwargs:
            tr_cur[:,:,:starti+i] = traj[:,:,:starti+i];
            x = np.hstack(tr_cur[0,:,:]);
            y = np.hstack(tr_cur[1,:,:]);
            z = np.hstack(tr_cur[2,:,:]);
            tracks.mlab_source.set(x=x,y=y,z=z);
        d['fig'].scene.disable_render=False;
        print("saving {}".format(names['out']));
        mlab.savefig(names['out'],size=(1280,1024));
    pass

def plot_single(inname,vlim,angle,
                **kwargs):
    print("loading mlab");
    import mayavi.mlab as mlab;
    initial_plot(mlab,inname,vlim,angle,
                 **kwargs);
    print('plotting');
    mlab.show();

if __name__=="__main__":
    #reading in arguments
    opts=docopt(__doc__,help=True);
    vmin = float(opts['--min']) if opts['--min'] else -1;
    vmax = float(opts['--max']) if opts['--max'] else 23.5;
    kwargs = {};
    kwargs['zeros']=[conv(opts['--zero-x'],func=eval,default=(None,None)),
                     conv(opts['--zero-y'],func=eval,default=(None,None)),
                     conv(opts['--zero-z'],func=eval,default=(None,None))];
    kwargs['zeros'] = None if kwargs['zeros'] == [(None,None)]*3 else kwargs['zeros'];
    angle = [ float(opts['--polar']),
              float(opts['--azimuth']),
              float(opts['--roll']) ];
    kwargs['clabel'] = opts['--clabel'];
    traj = opts['--trajectories'];
    if  opts['--cutplane']:
        kwargs['cutplane'] = int(opts['--cutplane']);
    if not opts['--plot-single']:
        if '*' not in opts['IN_FORMAT'] or '*' not in opts['OUT_FORMAT']:
            print('point format should contain \'*\'');
            print(usage);
            quit(1);
        lownum=int(opts['<lownum>']);
        highnum=int(opts['<highnum>']);
        numrange = range(lownum,highnum+1);
        if opts['<step>']:
            numrange = numrange[::int(opts['<step>'])];
        fmt='{{:0>{}}}'.format(len(opts['<highnum>']))
        out_fmt=re.sub(r'\*',fmt,opts['OUT_FORMAT']);
        in_fmt =re.sub(r'\*',fmt,opts['IN_FORMAT']);
        label_fmt=re.sub(r'\*',fmt,opts['LABEL_FORMAT']);
        outnames = [out_fmt.format(i) for i in numrange];#this so I don't need to check for None
        innames  = [in_fmt.format(i) for i in numrange];
        labels   = [label_fmt.format(i) for i in numrange];
        stuff = {'out':outnames,'in':innames,'label':labels};
        if traj:
            if type(traj) is str:
                if '*' not in traj:
                    print('point format should contain \'*\'');
                    print(usage)
                    quit(1);
            else:
                points_fmt = 'points*.pt'
            points_fmt=re.sub(r'\*',fmt,traj);
            points_files = (points_fmt.format(numrange[0]),
                            points_fmt.format(numrange[-1]));
            kwargs.update({'traj-files':points_files});
            print kwargs;
            # points = [points_fmt.format(i) for i in numrange];
            # stuff.update({'traj':points});
        files = [ dict(zip(stuff.keys(),d)) for d in zip(*stuff.values())];
        plot(files,(vmin,vmax),angle,**kwargs);
    else:
        plot_single(opts['INFILE'],(vmin,vmax),angle,**kwargs);
    pass;
pass;
