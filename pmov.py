#!/usr/bin/env python2
'''
Read a pmovie. In the absense of an output name, output
the name based on the frame step.

Usage:
    pmov.py [options] <input> [<output>]
    pmov.py [options] [--hdf | -H ] <input> <output>

Options:
    --help -h      Output this help.
    --sort -s      Sort the pmovies by IPP.
    --hdf -H       Output to hdf5 instead of to a pickle file.
                   The group will be based on the step.
    --zip -z       Use compression for hdf5.
    --verbose -v   Be verbose.
    --lock=L -l L  Specify a lock file for synchronized output for hdf5.
'''
import lspreader as rd;
import h5py as h5;
from misc import dump_pickle, h5w, mkvprint;
from docopt import docopt;
import numpy as np;
import fasteners;

opts = docopt(__doc__,help=True);
vprint = mkvprint(opts);


def sortframe(frame):
    '''
    sorts particles for a frame
    '''
    d = frame['data'];
    sortedargs = np.lexsort([d['xi'],d['yi'],d['zi']])
    d = d[sortedargs];
    frame['data']=d;
    return frame;

def hdfoutput(outname, frames, dozip=False):
    '''Outputs the frames to an hdf file.'''
    with h5.File(outname,'a') as f:
        for frame in frames:
            group=str(frame['step']);
            h5w(f, frame, group=group,
                compression='lzf' if dozip else None);
frames=rd.read(opts['<input>']);
if opts['--sort']:
    vprint("sorting...");
    frames[:] = [sortframe(frame) for frame in frames];
    vprint("done");
if opts['--hdf']:
    output = lambda :hdfoutput(opts['<output>'], frames, opts['--zip']);
    if opts['--lock']:
        output = fasteners.interprocess_locked(opts['--lock'])(output);
    output();
else:
    if not opts['<output>']:
        for frame in frames:
            outname = "{}.{}".format(opts['<input>'],frame['step']);
            dump_pickle(outname, frame);
    else:
        dump_pickle(opts['<output>'], frames);
