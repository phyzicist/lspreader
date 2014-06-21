#!/usr/bin/env python2

import lspreader as rd;
import cPickle;
import sys;

def main():
    usage="usage: python2 scalar.py <input> <output>"
    if len(sys.argv) != 3:
        print(usage);
        exit();
    with rd.LspOutput(sys.argv[1]) as f:
        print(f.header);
        print("reading in data");
        data = f.get_data(pool_size=24,lazy=False);
    print("selecting data");
    print("dumping");
    with open(sys.argv[2],"wb") as f:
        cPickle.dump(data,f,2);
    pass;

if __name__ == '__main__':
    main();
