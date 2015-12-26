#!/usr/bin/env python2
'''
An attempt to formalize the parsing of pmovies.

Usage: ./pmovie-prep.py [--help | -h]
'''
import subprocess;
import numpy as np;
import re;
from docopt import docopt;
docopt(__doc__);
def call(cmd):
    return subprocess.check_output(cmd).split('\n');
ls = call(('ls',))
# filtering for pmovies
pmovierx=re.compile(r"pmovie([0-9]+).p4$");
ls = [(s,int(pmovierx.match(s).group(1))) for s in ls if pmovierx.match(s)];
ls.sort(key=lambda i: i[1]);
if len(ls) == 0:
    print("I see no pmovies, exiting.");
    exit(0);
filesstr = [i[0] for i in ls];
#obtaining sizes
du = call(['du','-b']+filesstr);

sizerx = re.compile(r"^([0-9]+).*");
du = [int(sizerx.match(s).group(1)) for s in du if sizerx.match(s)];
#now, first, we divide into continguous groups;
sizes = np.array(du);
ds = sizes[1:] - sizes[:-1];
breaks  = np.where(ds != 0)[0];
if len(breaks) == 0:
    exit();
breaks += 1;#gettings the indices of breaks.
for i in breaks:
    print('{},{}'.format(ls[i-1][0],ls[i][0]));


