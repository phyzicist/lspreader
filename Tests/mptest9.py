# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 10:50:47 2016

Modifed from https://docs.python.org/2/library/multiprocessing.html

Demonstrates how you can assign tasks and then gather the results

@author: Scott
"""

from multiprocessing import Pool, TimeoutError
import time
import os

def f(x):
    return x*x

def f2(x, y):
    time.sleep(2)
    return x*y, x+y

if __name__ == '__main__':
    pool = Pool(processes=4)              # start 4 worker processes

    # evaluate "f(20)" asynchronously
    print "Assigning single task"
    res = pool.apply_async(f2, (3,4))      # runs in *only* one process
    print "Waiting for single response"
    print res.get(timeout=3)              # prints "12"

    args = [(3,4), (2,2), (1,10), (12,2), (13,1)]
    nargs = len(args)
    reslist = [None]*nargs
    print "Assigning multiple tasks"
    for i in range(nargs):
        reslist[i] = pool.apply_async(f2, args[i])      # runs in *only* one process
    
    print "Waiting for multiple replies"
    for i in range(nargs):
        print reslist[i].get(timeout=3)              # prints "12"
