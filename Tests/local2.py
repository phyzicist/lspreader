from h5stitch2D import *
from lspreader2 import *
import os
import gzip

#folder = r'/home/feister.7/lsp/runs/run2_lastest/'
folder = r'/tmp/ngirmang.1-2d-nosolid-151113/'
#fns = getfns(folder, ext='p4.gz')
#print len(getfns(folder, ext='p4.gz'))
path1 = os.path.join(folder,'flds6.p4.gz')
path2 = os.path.join(folder, 'flds8.p4.gz')

with gzip.open(path1, 'rb') as f:
	header1 = get_header(f)
	time1 = header1['timestamp']

with gzip.open(path2, 'rb') as f:
	time2 = get_header(f)['timestamp']

dt = time2 - time1 # timestep in ns
dt_secs = dt * 10**(-9) # timestep in seconds
dt_fs = dt * 10**6
print dt_fs

c = 3e8 # Speed of light in m/s
wl = 0.8e-6 # Wavelength of laser in m
t_cycle = wl/c # time duration of a laser cycle, in seconds

print t_cycle * 10**15

dt_cycles = dt_secs / t_cycle # Timestep in laser cycles
print dt_cycles

#h5path = r'/home/feister.7/lsp/runs/greg_run/fields2D.hdf5'
#h5path = h5fields2Dser(folder, h5path=h5path)

