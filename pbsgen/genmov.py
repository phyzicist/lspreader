#!/usr/bin/env python
'''
Make a mov PBS script.

Usage:
  mymov.py [options] <start> <stop> [<step>]
  mymov.py [options] (--stdin|-I)

Options:
  --hours=HOURS -H HOURS          Set the number of hours. [default: 02]
  --minutes=MINUTES -m MINUTES    Set the number of hours. [default: 00]
  --server=SERVER -n SERVER       For this server. [default: oakley]
  --post=POST -p POST             Set the postfix for the name of the job.
  --subdivisions=S -s S           Set subdivisions. [default: 1]
  --stdin -I                      Read numbers in from stdin.
  --clear                         Cut out a section.
  --indir=INDIR -i INDIR          Files are prefixed by INDIR [default: .]
  --outdir=OUTDIR -o OUTDIR       Files are prefixed by OUTDIR [default: .]
  --labeltype=T                   Set the label type [default: $I]
'''
from docopt import docopt;
import sys;
opts=docopt(__doc__,help=True);

header='''#PBS -l walltime={hours}:{min}:00
#PBS -l {server}
#PBS -S /nfs/04/osu8362/bin/bash
#PBS -j oe
#PBS -N mov{post}

source $HOME/.bashrc
source $HOME/conda

PBSBASE=mov{post}

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE | uniq >>${{PBSBASE}}.log


module load virtualgl
DNUM=`vncserver 2>&1 | grep : | head -1 | sed 's/.*\\(:.\\).*/\\1/'`
sleep 10
export DISPLAY=$DNUM
echo "opened $DISPLAY">>${{PBSBASE}}.log
'''
loop='''for ((i={first};i<={second};i+={step})); do
    while [ `pgrep -u $USER python2 | wc -l ` -gt 7 ]; do
       echo waiting...>>${{PBSBASE}}.log;
       sleep 5;
    done;
    I=`printf "%0{padlen}d" $i`
    L={label}
    echo RhoN10_${{I}}.H>>${{PBSBASE}}.log
    ./render3d.py -s "{indir}/RhoN10_${{I}}.H" "{outdir}/${{I}}.png" "t = ${{L}} fs" \
    {clear} --new-ctf &>>${{PBSBASE}}.log &
done;
while [ `pgrep -u $USER python2 | wc -l ` -gt 0 ]; do
    echo waiting on exiting ones>>${{PBSBASE}}.log;
    sleep 25;
done;
vncserver -kill $DNUM
killall bonobo-activation-server
'''
clear="--zero-x='(0,50)' --zero-y='(50,100)' --zero-z='(0,50)' --otf=3"
single="./sclr.py -xyz sclr{num}.p4 RhoN10 RhoN10_{num}.H>>${{PBSBASE}}.log";

label = opts['--labeltype']

if opts['--stdin']:
    r = [single.format(num=line.strip()) for line in sys.stdin.readlines()];
    out = header.format(hours=opts['--hours'],min=opts['--minutes'],ppn=opts['--ppn'],post='')+'\n'.join(r);
    print(out);
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
    post = name_pad.format(i);
    name = ('mov{}.pbs').format(post);
    curout = out.format(hours=opts['--hours'],min=opts['--minutes'],
            server="nodes=1:ppn=12:gpus=2:vis" if opts['--server'] == "oakley" else "nodes=1:ppn=8:gpu",
            first=f,second=s,step=opts['<step>'],
            label=label,
            post=post,clear = clear if opts['--clear'] else "",
            padlen=padlen,indir=opts['--indir'],outdir=opts['--outdir'],);
    with open(name,'w') as f:
        f.write(curout);
    pass;
pass;
