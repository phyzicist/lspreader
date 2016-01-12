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
    --outdir=OUT               Specify an output directory for pbs scripts [default: .]
    --hdf                      Use hdf5 as output.
    --ramses-node=NODE         Submit to a specific ramses node.
    --filelist=LS -l LS        Supply a directory listing if the directory is
                               unreachable for this script.
    --dotlsp=DOT -. DOT        Supply a .lsp file if it is unreachable for this
                               script.
'''
import subprocess;
import numpy as np;
import re;
from docopt import docopt;

opts = docopt(__doc__,help=True);
workdir = opts['<workdir>'];

def chunks(l,n):
    #http://stackoverflow.com/a/3226719
    #...not that this is hard to understand.
    return [l[x:x+n] for x in xrange(0, len(l), n)];

def call(cmd):
    return subprocess.check_output(cmd).split('\n');

def filelines(fname,strip=False):
    with open(fname,'r') as f:
        lines = f.readlines();
    if strip:
        lines[:] = [line.strip() for line in lines]
    return lines;

if not opts['--filelist']:
    ls = call(('ls',workdir))
else:
    ls = filelines(opts['--filelist'],strip=True);
# filtering for pmovies
pmovierx=re.compile(r"pmovie([0-9]+).p4$");
lspmovies = [(s, int(pmovierx.match(s).group(1)))
             for s in ls if pmovierx.match(s)];
lspmovies.sort(key=lambda i: i[1]);
if len(ls) == 0:
    raise ValueError("I see no pmovies in {}.".format(workdir));
filesstr = [i[0] for i in lspmovies];
#reading .lsp file to find information
if not opts['--dotlsp']:
    dotlsprx = re.compile(r".*\.lsp");
    dotlsp = [s for s in ls if dotlsprx.match(s)][0];
else:
    dotlsp = opts['--dotlsp'];
dotlsp = filelines(dotlsp,strip=True);

#Figure out the dimensionality of the simulation.
#Another possiblity is to read the output dir which
#has compiler flags. This will work for now as long
#as the .lsp file is cleanish. Basically search for the
#Grid and look for the *-cells entries
cellrxs=[re.compile(r"^{}-cells +[0-9]+".format(dim))
         for dim in ['x','y','z']]
dims =  [ np.array([rx.match(line) for line in dotlsp]).max()
          for rx in cellrxs ]
dims = [dim for dim,isdim in zip(['x','y','z'],dims) if isdim];
dims_flag = ' -{} '.format(''.join(dims));

#providing for different settings which varies
#by server.
server_settings = {
    #          maxproc ppn  hours  mins
    'ramses': [23,     48,   999,      0],
    'glenn':  [8,       8,    48,      0]
};
try:
    maxproc,ppn,hours,mins = server_settings[opts['--server']];
except KeyError as k:
    print('Invalid server "{}"'.format(k));
    exit(1);
#determining subdivisions
nodes = int(opts['--nodes']);
filespernode = len(filesstr)/nodes;
if len(filesstr)%nodes > 0: filespernode+=1;
pbses = chunks(filesstr, filespernode);
postfmt = '{{0:{}}}'.format(len(str(nodes)));

#template for PBS script output
template='''
#PBS -l walltime={hours}:{mins}:00
#PBS -l nodes=1:ppn={ppn}{ramsesnode}
#PBS -S /bin/bash
#PBS -j oe
#PBS -N pmovie-conv-{post}

source $HOME/.bashrc
source $HOME/conda
LOGFILE=$PBS_O_WORKDIR/pmovie-conv-{post}.log
MAXPROC={maxproc}
cd {workdir}
>$LOGFILE

if [ ! -d pmovie-conv ]; then
    mkdir pmovie-conv;
fi

echo "processing first file...">>$LOGFILE
{firstfile}
for i in {filelist}; do
    while [ $(pgrep -f pmov.py  |  wc -l ) -gt $MAXPROC ]; do sleep 10; done; 
    echo "running $i">>$LOGFILE
    sleep 1;
    ./pmov.py {opts} -D pmovie-conv  --exp-d=./hash.d $i {outfile}&>>$LOGFILE&
done

while [ $(pgrep -f pmov.py | wc -l) -gt 0 ]; do
    echo "waiting for $(pgrep -f pmov.py | wc -l) process(es)">>$LOGFILE
    echo "call deq to end">>$LOGFILE
    sleep 10;
done
'''

firstfiletmpl="./pmov.py {opts} -D pmovie-conv --exp-first=./hash.d {firstfile} {outfile}&>>$LOGFILE"
for i,f in enumerate(pbses):
    post = postfmt.format(i);
    xopts=''+dims_flag; #for copy so we don't write to dims_flag
    #handle hdf output
    if opts['--hdf']:
        xopts+=' -H '
        outfile='pmovs.h5'
    else:
        outfile='';
    #ramses hack
    if opts['--server'] and opts['--ramses-node']:
        ramsesnode=':'+opts['--ramses-node'];
    else:
        ramsesnode='';
    if i==0:
        firstfile = firstfiletmpl.format(
            opts=xopts,
            firstfile=f[0],
            outfile=outfile);
        files=' '.join(f[1:]);
    else:
        firstfile='';
        files = ' '.join(f);
    out = template.format(
        hours=hours,
        mins=mins,
        ppn=ppn,
        ramsesnode=ramsesnode,
        post=post,
        workdir=workdir,
        firstfile=firstfile,
        filelist=files,
        opts=xopts+opts['--extra-opts'],
        maxproc=maxproc,
        outfile=outfile);
    with open(opts['--outdir']+'/pmovie-conv-'+post+'.pbs','w') as f:
        f.write(out);
