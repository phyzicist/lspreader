#!/usr/bin/env python
import numpy as np;
import cPickle;
import sys;

def read_data(filename):
    print('loading file {}'.format(filename));
    with open(filename,'rb') as f:
        d=cPickle.load(f);
    return d;

def plot(S,outname):
    print("loading mlab");
    import mayavi.mlab as mlab;
    mlab.options.offscreen = True;
    print("plotting");
    src=mlab.pipeline.scalar_field(S);
    vmin=20;vmax=24;
    v=mlab.pipeline.volume(src,vmin=vmin,vmax=vmax);    
    mlab.view(elevation=150,focalpoint='auto',distance='auto');
    mlab.scalarbar(object=v,title="log10 of number density");
    #Unfortunately, the volume renderer doesn't work out of the box.
    #The main issue is the colorbar. So, we do this instead.
    e=mlab.get_engine();
    module_manager = e.scenes[0].children[0].children[0];
    module_manager.scalar_lut_manager.use_default_range = False;
    module_manager.scalar_lut_manager.data_range = np.array([vmin, vmax]);
    mlab.text(0.075,0.875,name);
    mlab.savefig(outname);

if __name__=="__main__":
    usage="usage: ./plot3d.py <input> <output.png>";
    if len(sys.argv) != 3:
        print(usage);
        exit(1);
    global name;
    inname=sys.argv[1];
    name=sys.argv[1][:-2];
    outname=sys.argv[2];
    S=read_data(inname);
    plot(S,outname);
    print('plotted');
pass;
