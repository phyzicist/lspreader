# Run this first to create the hdf5 file
import h5py

with h5py.File('parallel_test.hdf5', 'w') as f:
    dset = f.create_dataset('test', (4,), dtype='i')
