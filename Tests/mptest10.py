# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 16:39:44 2016

Nonblocking recieve
Copied from: http://stackoverflow.com/questions/25001912/non-blocking-way-of-receiving-a-pickled-item-in-mpi4py

@author: Scott
"""

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
   data = {'a': 7, 'b': 3.14}
   time.sleep(6)
   comm.send(data, dest=1, tag=11)

elif rank == 1:
   while not comm.Iprobe(source=0, tag=11):
        print 'rank 1 Doing some work...'
        time.sleep(1)
   rdata = comm.recv(source=0, tag=11)
   print 'rank 1: got ', rdata
   rdata = comm.recv(source=2, tag=0)
   print 'rank 1: got ', rdata

elif rank == 2:
   data = {'a': 7, 'b': 3.24}
   time.sleep(3)
   comm.send(data, dest=1, tag=0)