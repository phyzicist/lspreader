'''
Search a file for the given indices.

Usage:
    ./search.py [options] <input> <indexfile> <output>

Options:
    --help -h                 Print this help.
'''
from docopt import docopt;
opts=docopt(__doc__,help=True);
from misc import readfile;
import numpy as np
indices=np.load(opts['<indexfile>']);
current=readfile(opts['<input>'],dumpfull);
found  =np.in1d(current['data']['hash'],indices);
#boolean array to indices
np.arange(len(indices));
