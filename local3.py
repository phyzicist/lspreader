import h5stitch2D as hs
from freqanalysis import *
import lspreader as rd
import numpy as np
import os

h5path = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir/test/fields2D.hdf5'

folder = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir'
outdir = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir/test/plots'
fns = hs.getfnsp4(folder)
print len(fns)


#for fld_id in ['Ex', 'Ez', 'By']:
#    data2 = freqanalyze(fns, fld_id = fld_id)

#h5path = hs.h5fields2D(folder)
#data = fa.loadData2(h5path)
#for fld_id in ['Ex', 'Ez', 'By']:
#    data2 = freqanalyze(fns, fld_id = fld_id, data = data)

#h5path = hs.h5fields2D(folder)
#for fld_id in ['Ex', 'Ez', 'By']:
#    data2 = freqanalyze(fns, fld_id=fld_id, datpath = h5path)

def subdir(folder, name):
    """ Make a subdirectory in the specified folder, if it doesn't already exist"""
    subpath = os.path.join(folder,name)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    return subpath

fld_ids = ['Ex', 'Ez', 'By']
data = hs.fields2D(fns, fld_ids = fld_ids)
for fld_id in fld_ids:

    data2 = freqanalyze(fns, fld_id = fld_id, data = data)

    # Give the subdirectory for this dataset a name, and make the directory.
    meantime = np.mean(data2['times_fs']) # Mean time, in fs
    label = "{:.0f}".format(round(meantime*10)).zfill(5) # Make the folder label be '00512' for t = 51.2342 fs
    
    flddir = subdir(outdir, fld_id) # E.g. make a directory called "[outdir]/Ex" if it doesn't exist
    timedir = subdir(flddir, label) # E.g. make a directory called "[outdir]/Ex/0127" if it doesn't exist, for time 12.7 fs

    pltdict = getPltDict(data2)
    h5path = saveData2(data2, folder=timedir)
    plotme(data2, folder=timedir, pltdict=pltdict)
