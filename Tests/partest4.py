# Copied from http://docs.h5py.org/en/latest/mpi.html
# Run via: mpiexec -n 4 python partest4.py
# Then, check how it did with: h5dump parallel_test.hdf5
from mpi4py import MPI
import h5py

rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)

f = h5py.File('parallel_test.hdf5', 'w', driver='mpio', comm=MPI.COMM_WORLD)

dset = f.create_dataset('test', (4,), dtype='i')
dset[rank] = rank

f.close()
