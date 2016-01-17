# Copied from http://docs.h5py.org/en/latest/mpi.html
# Run via: mpiexec -n 4 python partest5.py

from h5stitch2D import h5fields2Dpar

folder = r'/media/sf_temp/lastest/run2_lastest longer/testdir'

h5path = h5fields2Dpar(folder) # Does not do any compression.
