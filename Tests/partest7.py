# Based on example from: https://www.howtoforge.com/tutorial/distributed-parallel-programming-python-mpi4py/#-introduction
# Shows the right way to do parallel operations

from mpi4py import MPI

def chunkIt(seq, num):
    """Divide list 'seq' into 'num' (nearly) equal-sized chunks. Returns a list of lists; some empty lists if num is larger than len(seq)."""
    # Copied directly from: http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

nprocs = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

comm = MPI.COMM_WORLD

if rank == 0:
   full_list = ['a','b','c','d','e','f','g','h','i','j'] # Define the full list
   data = chunkIt(full_list, nprocs)
   print 'scattering data',data
else:
   data = None
data = comm.scatter(data,root=0)
print 'rank',rank,'has data: ', data

for i in range(len(data)):
    data[i] = data[i] + str(rank)
 
new_data = comm.gather(data, root=0)
if rank == 0:
   print  'master collected: ', new_data
