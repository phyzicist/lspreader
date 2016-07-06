# Assuming LSP option forcing one movie frame per .p4 file

import numpy as np

def firsthash(frame,dims,removedupes=False):
    '''
	[Copied from Gregory Ngirmang's lspreader.
	https://github.com/noobermin/lspreader]
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
	[Copied from Gregory Ngirmang's lspreader.
	https://github.com/noobermin/lspreader]
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


frames = rd.read_movie2('pmovie1.p4')
frame = frames[0]
hashd = firsthash(frame, ['xi','zi','yi'], removedupes=True)
hashes = genhash(frame, hashd, removedupes=True)
ix_sort = np.argsort(hashes)
hash_ref = hashes[ix_sort] # Sorted list of hashes
data_sort = frame['data'][ix_sort]

for fn in ['pmovie20.p4','pmovie35.p4']:
	frames = rd.read_movie2(fn)
	frame = frames[0]
	frame = addhash(frame, hashd, removedupes=True)
	framedat_sort = np.sort(frame['data'], order='hash')
	# Iterate over the framereference; if it matches the hash, include it. Otherwise, Call it bad.
	
	