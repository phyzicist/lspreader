# Run this second to add all the files to the database
from mpi4py import MPI
import h5py
from time import sleep

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

try:
    with h5py.File('parallel_test.hdf5', 'r+') as f:
        print "Successfully opened file on processor " + str(rank)
        f['test'][rank] = rank + 10
        sleep(3)
except:
    print "Failed to open file on processor " + str(rank)
