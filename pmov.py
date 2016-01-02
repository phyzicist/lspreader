#!/usr/bin/env python2
'''
Read a pmovie. In the absense of an output name, output
the name based on the frame step.

Usage:
    pmov.py [options] <input> [<output>]
    pmov.py [options] --hdf <input> <output>

Options:
    --help -h      Output this help.
    --sort -s      Sort the pmovies by IPP. Implies --unique.
    --unique -u    Take unique particles.
    --hdf -H       Output to hdf5 instead of to a pickle file.
                   The group will be based on the step.
    --zip -z       Use compression for hdf5.
'''
import lspreader2 as rd;
from misc import dump_pickle, h5w;
from docopt import docopt;

opts = docopt(__doc__,help=True);
frames=rd.read(opts['<input>']);
if opts['--sort']: opts['--unique']=True;

for i,frame in enumerate(frames):
    d = frame['data'];
    if opts['--unique']:
        _,uniq = np.unique(d[['xi','yi','zi']]);
        frame = frame[uniq];
    if opts['--sort']:
        sortedargs = np.lexsort([d['xi'],d['yi'],d['zi']])
        frame = frame[sortedargs];
    frame['data']=d;
    frames[i]=frame;

if opts['--hdf']:
    if opts['--zip']:
        c = 'lzf';
    else:
        c = None;
    with h5.File(opts['<output>'],'a') as f:
        for frame in frames:
            group=str(frame['step']);
            h5w(f, frame,
                group=group, compression=c);
else:
    if not opts['<output>']:
        for frame in frames:
            outname = "{}.{}".format(opts['<input>'],frame['step']);
            dump_pickle(outname, frame);
    else:
        dump_pickle(opts['<output>'], frames);

