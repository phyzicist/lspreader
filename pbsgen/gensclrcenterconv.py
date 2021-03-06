#!/usr/bin/env python
'''
Make a sclrcenterconv PBS script.

Usage:
  gensclrcenterconv.py [options] <start> <stop> [<step>]
  gensclrcenterconv.py [options] (--stdin|-I)

Options:
  --hours=HOURS -H HOURS            Set the number of hours. [default: 02]
  --minutes=MINUTES -m MINUTES      Set the number of hours. [default: 00]
  --ppn=PPN -n PPN                  Set the ppn. [default: 1]
  --post=POST -p POST               Set the postfix for the name of the job. [default: ]
  --subdivisions=S -s S             Set subdivisions. [default: 1]
  --stdin -I                        Read numbers in from stdin.
  --var=VAR -v VAR                  Read the following variables in a list. [default: ["RhoN10"] ]
  --outdir=OUTDIR -o OUTDIR         Output files to directory OUTDIR [default: .]
  --extra-options=OPTS -x OPTS      Send extra options to every command. [default: ]
  --server=SERVER                   Make for server. [default: glenn]
'''
from docopt import docopt;
import sys;
opts=docopt(__doc__,help=True);

header='''#PBS -l walltime={hours}:{min}:00
#PBS -l nodes=1:ppn={ppn}{extra_lopts}
#PBS -S /nfs/04/osu8362/bin/bash
#PBS -j oe
#PBS -N sclrcenterconv{post}

source $HOME/.bashrc
source $HOME/conda

PBSBASE=sclrcenterconv{post}

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE | uniq >${{PBSBASE}}.log
'''
loop='''for ((i={first};i<={second};i+={step})); do
    I=`printf "%0{padlen}d" $i`
    echo sclr${{I}}.p4&>>${{PBSBASE}}.log
    ./sclrcenter.py --X --half sclr${{I}}.p4 {varpairs} {opts} &>>${{PBSBASE}}.log
done;
'''
single="./sclrcenter.py --X --half sclr{num}.p4 {varpairs} {opts} &>>${{PBSBASE}}.log";

varpairs=' '.join(['{var} {outdir}/{var}-{{I}}.s'.format(var=s,outdir=opts['--outdir'])
                   for s in  eval(opts['--var'])]);

if opts['--server'] == 'ramses':
    extra_lopts=':hedp';
else:
    extra_lopts='';
    
if opts['--stdin']:
    r = [single.format(num=line.strip(),
                       varpairs=varpairs.format(I=line.strip()),
                       opts=opts['--extra-options'])
                       for line in sys.stdin.readlines()];
    body = '\n'.join(r);
    head = header.format(hours=opts['--hours'],min=opts['--minutes'],ppn=opts['--ppn'],
                        post=opts['--post'],extra_lopts=extra_lopts);
    print(head+body);
    quit();


first = int(opts['<start>']);
stop = int(opts['<stop>']);
step = int(opts['<step>']) if opts['<step>'] else 1;
r = range(first,stop+step,step);
subdiv = int(opts['--subdivisions']);
name_pad = '{{:0{}}}'.format(len(opts['--subdivisions']));
subdivs = zip(r[::len(r)/subdiv],r[len(r)/subdiv-1::len(r)/subdiv][:-1]+[r[-1]]);
padlen = len(str(r[-1]));
out = header+loop;
for i,cur in enumerate(subdivs):
    f,s=cur;
    post = opts['--post']+name_pad.format(i);
    name = ('sclrcenterconv{}.pbs').format(post);
    curout = out.format(hours=opts['--hours'],min=opts['--minutes'],
            ppn=opts['--ppn'],first=f,second=s,step=step,
            varpairs=varpairs.format(I="${I}"),
            post=post, padlen=padlen,
            opts=opts['--extra-options'],
            extra_lopts=extra_lopts);
    with open(name,'w') as f:
        f.write(curout);
    pass;
pass;
