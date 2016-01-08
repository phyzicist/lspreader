from h5stitch2D import *
#folder = r'/home/feister.7/lsp/runs/run2_lastest/'
folder = r'/tmp/ngirmang.1-2d-nosolid-151113/'
print len(getfns(folder, ext='p4.gz'))
h5path = r'/home/feister.7/lsp/runs/greg_run/fields2D.hdf5'
h5path = h5fields2Dser(folder, h5path=h5path)

