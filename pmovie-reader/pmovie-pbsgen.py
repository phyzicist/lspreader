#!/usr/bin/env python2
'''
An attempt to formalize the parsing of pmovies.

Usage: 
    ./pmovie-prep.py [options] <workdir>

Options:
    --help -h                  Print this help.
    --nodes=N -n N             Number of utilized nodes [default: 1]
    --server=SERVER            Tune for this server [default: ramses]
    --extra-opts=OPTS -x OPTS  Add these opts to pmov.py [default: ]
    --outdir=OUT               Specify an output directory [default: .]
'''
import subprocess;
import numpy as np;
import re;
from docopt import docopt;
from misc import chunks;
from pprint import pprint;

opts = docopt(__doc__,help=True);

def call(cmd):
    return subprocess.check_output(cmd).split('\n');
ls = call(('ls',))
# filtering for pmovies
pmovierx=re.compile(r"pmovie([0-9]+).p4$");
lspmovies = [(s,int(pmovierx.match(s).group(1))) for s in ls if pmovierx.match(s)];
drx = re.compie(r"pmovie([0-9]+).p4.d$");
lsds = [int(drx.match(s).group(1)) for s in ls if drx.match(s)];
#O(N^2), although N is only ~2000, we can do better.
lspmovies = [i for i in lspmovies if i[1] not in lsds];
lspmovies.sort(key=lambda i: i[1]);
if len(ls) == 0:
    print("I see no pmovies, exiting.");
    exit(1);
filesstr = [i[0] for i in lspmovies];

#providing for different settings which varies
#by server.
server_settings = {
    #           ppn  hours  mins
    'ramses': [48,   999,      0],
    'glenn':  [ 8,    48,      0]
};
try:
    ppn,hours,mins = server_settings[opts['--server']];
except KeyError as k:
    print('Invalid server "{}"'.format(k));
    exit(1);
#determining subdivisions
nodes = int(opts['--nodes']);
filespernode = len(filesstr)/nodes;
if len(filesstr)%nodes > 0: filespernode+=1;
pbses = chunks(filesstr, filespernode);
#subdivs = nodes*ppn;
#filesperproc = len(filesstr)/subdivs;
#if len(filesstr)%subdivs > 0:  filesperproc += 1;
#pbses=chunks(filesstr, filesperproc);

workdir = opts['<workdir>'];
postfmt = '{{0:{}}}'.format(len(str(nodes)));

#this is the header
template='''
#PBS -l walltime={hours}:{mins}:00
#PBS -l nodes=1:ppn={ppn}
#PBS -S /bin/bash
#PBS -j oe
#PBS -N pmovie-conv-{post}

source $HOME/.bashrc
source $HOME/conda
LOGFILE=pmovie-conv-{post}.log
cd {workdir}

for i in "{filelist}"; do
    while [ $(pgrep pmov.py  |  wc -l ) -eq {ppn} ]; do sleep 10; done; 
    ./pmov.py -s {opts} $i {outdir}/$i.d &>>$LOGFILE&
done

while [ $(pgrep pmov.py | wc -l) -gt 0 ]; do
    echo "waiting for $(pgrep pmov.py | wc -l) process(es)">>$LOGFILE
    echo "call deq to end">>$LOGFILE
    sleep 10;
done
'''
pbsouts=[];

for i,f in enumerate(pbses):
    post = postfmt.format(i);
    files = ' '.join(f);
    out = template.format(
        ppn=ppn,
        hours=hours,
        mins=mins,
        post=post,
        workdir=workdir,
        filelist=files,
        opts=opts['--extra-opts'],
        outdir=opts['--outdir']);
    with open('pmovie-conv-'+post+'.pbs','w') as f:
        f.write(out);

        
