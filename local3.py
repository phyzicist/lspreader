import h5stitch2D as hs
from freqanalysis import *
import lspreader as rd
import numpy as np
import os
from multiprocessing import Pool

h5path = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir/test/fields2D.hdf5'




#for fld_id in ['Ex', 'Ez', 'By']:
#    data2 = freqanalyze(fns, fld_id = fld_id)

#h5path = hs.h5fields2D(folder)
#data = fa.loadData2(h5path)
#for fld_id in ['Ex', 'Ez', 'By']:
#    data2 = freqanalyze(fns, fld_id = fld_id, data = data)

#h5path = hs.h5fields2D(folder)
#for fld_id in ['Ex', 'Ez', 'By']:
#    data2 = freqanalyze(fns, fld_id=fld_id, datpath = h5path)


#fld_ids = ['Ex', 'Ez', 'By']
p4dir = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir/test'
divsp = 1
outdir = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir/test/plots'

freqFull(p4dir, outdir, nbatch = 20, npool = 4)
