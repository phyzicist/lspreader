#!/usr/bin/env python2
'''
Extract data directly.

Usage: read.py <input> <output>
'''
import lspreader2 as rd;
from misc import dump_pickle;
from docopt import docopt;

opts = docopt(__doc__,help=True);
dump_pickle(
    opts['<output>'],
    rd.read(opts['<input>']));

