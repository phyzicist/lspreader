# Copied from http://docs.h5py.org/en/latest/mpi.html
# Run via: mpiexec -n 4 python partest5.py

import h5py
import os
import shutil

folder = r'/media/sf_temp/lastest/run2_lastest longer/testdir'

h5path1 = os.path.join(folder,r'fields2D.hdf5')
h5path2 = os.path.join(folder,r'fields2D_comp.hdf5')

# Assume all keys are datasets, and copy them into the other hdf5
with h5py.File(h5path1, 'r') as f1:
    with h5py.File(h5path2, 'w') as f2:
        for k in f1.keys():
            f2.create_dataset(k, data=f1[k][...], compression='gzip', compression_opts=4)

shutil.move(h5path2, h5path1)
