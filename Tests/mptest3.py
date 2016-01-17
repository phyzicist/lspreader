from multiprocessing import Pool
import numpy as np

def myfun(a):
    return a, a+2

if __name__ == '__main__':
    p = Pool(4)
    fns = [1,1.01,1.02]
    c, d, e = p.map(myfun, fns)
    print "C:", c
    print "D:", d
    print "E:", e

