# Copied from http://docs.h5py.org/en/latest/mpi.html
# Run via: mpiexec -n 4 python partest4.py
# Then, check how it did with: h5dump parallel_test.hdf5
from mpi4py import MPI
import h5py
import numpy as np

nprocs = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

comm = MPI.COMM_WORLD

if rank==0:
    f = h5py.File('parallel_test.hdf5', 'w')
    dset = f.create_dataset('test', (4,), dtype='i')

testarr = np.ones(rank)*rank

full_list = comm.gather(rank, root=0)
flds = comm.gather(testarr, root=0) # Could speed this up, but do I need to?

if rank==0:
    print full_list
    for i in range(len(full_list)):
        dset[i]=full_list[i]
    f.close()
    print flds
