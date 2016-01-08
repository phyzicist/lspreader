#!/usr/bin/env python2
'''
Read a pmovie. In the absense of an output name, output
the name based on the frame step.

Usage:
    pmov.py [options] <input> [<output>]
    pmov.py [options] [--hdf | -H ] <input> <output>

Options:
    --help -h         Output this help.
    --sort -s         Sort the pmovies by IPP.
    --hdf -H          Output to hdf5 instead of to a pickle file.
                      The group will be based on the step.
    --zip -c          Use compression for hdf5.
    --npy -n          Use numpy.save to save.
    --verbose -v      Be verbose.
    --lock=L -l L     Specify a lock file for synchronized output for hdf5.
    --exp-d=DFILE     Experimental hashing. It might work only for uniform
                      grids. Specify specification file used to generate hash
                      as DFILE.
    --exp-first=DFILE Experimental hashing, see above. Use this file as the
                      first file to generate the hash specification from and
                      output to DFILE.
    --X -x            Use X as a spatial dimension. Similar options below are
                      for Y and Z. If none are passed, assume 3D cartesian.
    --Y -y            See above.
    --Z -z            See above.
'''
import lspreader as rd;
import h5py as h5;
from misc import dump_pickle, h5w, mkvprint,  readfile;
from docopt import docopt;
import numpy as np;
import numpy.lib.recfunctions as rfn;

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

def firsthash(frame,dims,removedupes=False):
    '''
    Hashes the first time step. Only will work as long as
    the hash can fit in a i8.

    Parameters:
    -----------
      frame : first frame.
      dims  :  iterable of strings for dimensions.

    Keywords:
    ---------
      removedups: specify duplicates for the given frame.
    
    Returns a dictionary of everything needed
    to generate hashes from the genhash function.
    
    '''
    #hashes must have i8 available
    #overwise, we'll have overflow
    def avgdiff(d):
        d=np.sort(d);
        d = d[1:] - d[:-1]
        return np.average(d[np.nonzero(d)]);
    ip    = np.array([frame['data'][l] for l in dims]).T;
    avgdiffs = np.array([avgdiff(a) for a in ip.T]);
    mins  = ip.min(axis=0);
    ips = (((ip - mins)/avgdiffs).round().astype('i8'))
    pws  = np.floor(np.log10(ips.max(axis=0))).astype('i8')+1
    pws = list(pws);
    pw = [0]+[ ipw+jpw for ipw,jpw in
               zip([0]+pws[:-1],pws[:-1]) ];
    pw = 10**np.array(pw);
    #the dictionary used for hashing
    d=dict(labels=dims, mins=mins, avgdiffs=avgdiffs, pw=pw);
    if removedupes:
        hashes = genhash(frame,d,removedupes=False);
        #consider if the negation of this is faster for genhash
        uni,counts = np.unique(hashes,return_counts=True);
        d.update({'dupes': uni[counts>1]})
    return d;

def genhash(frame,d,removedupes=False):
    '''
    Generate the hashes for the given frame for a specification
    given in the dictionary d returned from firsthash.

    Parameters:
    -----------
      frame :  frame to hash.
      d     :  hash specification generated from firsthash.

    Keywords:
    ---------
      removedups: put -1 in duplicates
    
    Returns an array of the shape of the frames with hashes.
    '''
    ip = np.array([frame['data'][l] for l in d['labels']]).T;
    scaled = ((ip - d['mins'])/d['avgdiffs']).round().astype('i8');
    hashes = (scaled*d['pw']).sum(axis=1);
    #marking duplicated particles
    if removedupes:
        dups = np.in1d(hashes,d['dupes'])
        hashes[dups] = -1
    return hashes;
def addhash(frame,d,removedupes=False):
    '''
    helper function to add hashes to the given frame
    given in the dictionary d returned from firsthash.

    Parameters:
    -----------
      frame :  frame to hash.
      d     :  hash specification generated from firsthash.

    Keywords:
    ---------
      removedups: put -1 in duplicates
    
    Returns frame with added hashes, although it will be added in
    place.
    '''
    hashes = genhash(frame,d,removedupes);
    frame['data'] = rfn.rec_append_fields(
        frame['data'],'hash',hashes);
    return frame;

#script start. This garbage to let people like scott use
#it without having to call it script wise.
if __name__=='__main__':
    opts = docopt(__doc__,help=True);
    vprint = mkvprint(opts);
    dims=[]
    if opts['--X']: dims.append('x');
    if opts['--Y']: dims.append('y');
    if opts['--Z']: dims.append('z');
    if len(dims)==0:
        dims=['x','y','z'];
    #reading in using the reader.
    frames=rd.read(opts['<input>']);
    
    if opts['--sort']:
        vprint("sorting...");
        frames[:] = [sortframe(frame) for frame in frames];
        vprint("done");
    #experimental hashing
    if opts['--exp-first']:
        d=firsthash(frames[0],dims);
        dump_pickle(opts['--exp-first'], d);
        frames[:] = [addhash(frame,d) for frame in frames];
    elif opts['--exp-d']:
        d = readfile(opts['--exp-d'],dumpfull=True);
        frames[:] = [addhash(frame,d) for frame in frames];
    #outputting.
    if opts['--hdf']:
        import fasteners;
        output = lambda :hdfoutput(opts['<output>'], frames, opts['--zip']);
        if opts['--lock']:
            output = fasteners.interprocess_locked(opts['--lock'])(output);
        output();
    elif opts['--npy']:
        if not opts['<output>']:
            for frame in frames:
                outname = "{}.{}".format(opts['<input>'],frame['step']);
                np.save(outname, frame['data']);
        else:
            np.save(opts['<output>'], [frame['data'] for frame in frames]);
    else:
        if not opts['<output>']:
            for frame in frames:
                outname = "{}.{}".format(opts['<input>'],frame['step']);
                dump_pickle(outname, frame);
        else:
            dump_pickle(opts['<output>'], frames);
