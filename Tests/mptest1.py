from multiprocessing import Pool

# Copied from https://docs.python.org/2/library/multiprocessing.html

def f(x):
    return x*x

if __name__ == '__main__':
    p = Pool(5)
    print(p.map(f, [1, 2, 3]))
