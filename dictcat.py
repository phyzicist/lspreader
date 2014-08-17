#!/usr/bin/env python2
'''
Concatenate dictionaries of numpy arrays, shape (1,) only.

Usage:
  dictcat.py [options] <input>...

Options:
  -o OUTPUT --output=OUTPUT         Output to OUTPUT instead of stdout.
'''
import cPickle;
import numpy as np;
from docopt import docopt;

def main():
    opts = docopt(__doc__,help=True);
    names = opts['<input>'];
    d=[]
    for name in names:
        with open(name) as f:
            d.append(cPickle.load(f))
    #making into lists
    data=d[0];
    for k in data:
        data[k] = list(data[k]);
    #stringing together
    for i in d[1:]:
        for k in data:
            data[k].extend(list(i[k]));
    pass;
    del d;
    #renumpying
    for k in data:
        data[k] = np.array(data[k]);
    if opts['--output']:
        with open(opts['--output'],"wb") as f:
            cPickle.dump(data,f,2);
        pass;
    else:
        print(cPickle.dumps(data,2));
    pass;
if __name__ == '__main__':
    main();
