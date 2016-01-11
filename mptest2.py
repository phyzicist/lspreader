from multiprocessing import Pool
import numpy as np

# Copied from https://docs.python.org/2/library/multiprocessing.html

def f(x):
    return x*x

def mymean(x):
    return np.mean(x)

def fftpar(Ei, npool = 10):
    p = Pool(npool)
    return np.swapaxes(np.array((p.map(np.fft.rfft, Ei.swapaxes(0,2)))),0,2)

def fftser(Ei):
    Eout = np.zeros((6, 10000, 200))
    for i in range(10000):
        for j in range(200):
            Eout[:,i,j] = np.fft.rfft(Ei[:,i,j])
    return Eout

if __name__ == '__main__':
    #myarr = np.array([[1,1,2],[2,3,3],[3,5,6],[2,3,4]])
    Ein = np.ones((10,10000,200))

    #print myarr

    Epar = fftpar(Ein)

    print Epar.shape
