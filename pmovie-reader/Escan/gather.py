'''
Gather found particles. Searches for "found*.npy" files in the current directory. Works on linux only.

Usage:
    ./gather.py [options] <output>

Options:
    --help -h                 Print this help.
'''
from docopt import docopt;
opts=docopt(__doc__,help=True);
import subprocess
def call(cmd):
    return subprocess.check_output(cmd).split('\n');
ls = call(('ls'))
foundrx = re.compile(r"found.*.npy");
files   = [file for file in ls
           if foundrx.match(file)]
arrays = [np.load(file) for file in files];
good   = np.concatenate(arrays);
good = np.unique(good);
np.save(opts['<output>'],good);
