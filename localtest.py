from h5stitch2D import *
folder = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir/test'
#folder = r'/home/feister.7/lsp/runs/run2_lastest/'
#folder = r'/tmp/ngirmang.1-2d-nosolid-151113/'
print len(getfnsp4(folder))
#h5path = r'/home/feister.7/lsp/runs/greg_run/fields2D.hdf5'
h5path = h5fields2D(folder)

