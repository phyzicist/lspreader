# Run this second to add all the files to the database. FAILS!!!!
from mpi4py import MPI
import h5py
import h5stitch2D as h52D

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

folder = r'/media/sf_temp/lastest/run2_lastest longer/testdir'
try:
    h5path, idx_part = h52D.h5fields2Db(folder, nprocs = size, rank = rank)
    print "Succeeded in adding " + str(len(idx_part)) + " files on processor " + str(rank)
except:
    print "Failed to add files on processor " + str(rank)
