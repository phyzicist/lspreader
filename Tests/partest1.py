# Demo of splitting a list across several stacks
# run via "mpiexec -n 4 python partest1.py"

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


size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()


full_list = ['a','b','c','d','e','f','g','h','i','j'] # Define the full list
chunk_list = chunkIt(full_list,size) # Break the list into N chunks, where N = number of processors
part_list = chunk_list[rank] # Assign each processor one chunk of the list

if rank < 1:
    print "Full list: ", full_list
#    print "This list will be split across " + str(size) + " processors as follows:"

print 'Proc ' + str(rank) + " of " + str(size) + " (on " + name + "): ", part_list

