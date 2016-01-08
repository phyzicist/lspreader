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
ls = call(('ls',workdir))
# filtering for pmovies
pmovierx=re.compile(r"pmovie([0-9]+).p4$");
lspmovies = [(s,int(pmovierx.match(s).group(1)))
             for s in ls if pmovierx.match(s)];
drx = re.compile(r"pmovie([0-9]+).p4.d$");
lsds = [int(drx.match(s).group(1))
        for s in ls if drx.match(s)];
#O(N^2), although N is only ~2000, we can do better.
lspmovies = [i for i in lspmovies if i[1] not in lsds];
lspmovies.sort(key=lambda i: i[1]);
if len(ls) == 0:
    print("I see no pmovies, exiting.");
    exit(1);
filesstr = [i[0] for i in lspmovies];
dotlsprx = re.compile(r".*\.lsp");
dotlsp = [s for s in ls if dotlsprx.match(s)][0];
#reading .lsp file to find information
with open(dotlsp,'r') as f:
    dotlsp = f.readlines();
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
    #           ppn  hours  mins
    'ramses': [44,   999,      0],
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
cd {workdir}
>$LOGFILE
echo "processing first file...{firstfile}">>$LOGFILE
{firstfile}
for i in {filelist}; do
    while [ $(pgrep -f pmov.py  |  wc -l ) -eq {ppn} ]; do sleep 10; done; 
    echo "running $i">>$LOGFILE
    ./pmov.py {opts} -n --exp-d=./hash.d $i {outfile}&>>$LOGFILE&
done

while [ $(pgrep -f pmov.py | wc -l) -gt 0 ]; do
    echo "waiting for $(pgrep -f pmov.py | wc -l) process(es)">>$LOGFILE
    echo "call deq to end">>$LOGFILE
    sleep 10;
done
'''

firstfiletmpl="./pmov.py -n --exp-first=./hash.d {opts} {firstfile} {outfile}&>>$LOGFILE"
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
        outfile=outfile);
    with open(opts['--outdir']+'/pmovie-conv-'+post+'.pbs','w') as f:
        f.write(out);
