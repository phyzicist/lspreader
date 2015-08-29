#!/usr/bin/env python2
'''
Make two .d consistent such that any missing ipp's are filled in with NaN's on the other
file. Assumes the inputs are sorted with pmov-sort.py

Usage:
  pmov-consistent.py [options] <first> <second> <output>

Options:
  --help -h             You know.
  --X=XLIM -x XLIM      Use X with limits. [default: (-30e-4, 5e-4)]
  --Y=YLIM -y YLIM      Use Y with limits. [default: None]
  --Z=ZLIM -z ZLIM      Use Z with limits. [default: (-20e-4,20e-4)]
  --pad -p              Pad the numbers with 0's.
  --hack -H             Only output to the second.
  --verbose -v          Be verbose.
'''

from docopt import docopt;
import cPickle as pickle;
from itertools import izip,imap;
import numpy as np;
import re;
import misc as m;

opts=docopt(__doc__,help=True);
def _print(s):
    print(s);
if not opts['--hack']: raise NotImplementedError("General case not implemented yet.");
vprint = _print if opts['--verbose'] else lambda s: None;
lims = [eval(opts[s]) for s in ["--X","--Y","--Z"]];
dims = [name for i,name in zip(lims,['xi','yi','zi']) if i is not None];
lims  = [lim for lim in lims if lim is not None];

if lims == [None,None,None]:
    raise NotImplementedError("Limitless simulations are not implemented");
vprint('reading first file {}'.format(opts['<first>']));
first =m.readfile(opts['<first>'],dumpfull=True);
vprint('reading second file {}'.format(opts['<second>']));
second=m.readfile(opts['<second>'],dumpfull=True);

data1 = first[-1]['data'];
data2 = second[0]['data'];
pi1   = np.array([data1[dim].astype("f4") for dim in dims]).T
pi2   = np.array([data2[dim].astype("f4") for dim in dims]).T
vprint("starting comparison");
vprint("size of first:{}\nsize of second:{}".format(len(pi1),len(pi2)));

def getfills1(pi1,pi2):
    #indices along the two arrays
    i,j=0,0;
    #nan insertions, rows for data and indices
    nanrows1=[]; nani1=[];
    nanrows2=[]; nani2=[];    
    while True:
        if i > len(pi1) or j > len(pi2): break;
        iend,jend = len(pi1),len(pi2);
        if iend - i > jend - j:
            iend -= (iend - i - jend + j);
        elif iend - i > jend - j:
            jend -= (jend - j - iend+ i);
        check = np.isclose(pi1[i:iend],pi2[j:jend]);
        addedcheck = check.T[0] & check.T[1];
        vprint("i,j={},{}".format(i,j));
        worst = np.argmin(addedcheck);
        if addedcheck[worst] == True:
            i=iend;
            j=jend;
            break;
        vprint("worst-i: {},{}".format(worst+i,check[worst]));
        #looking at the mismatch from the first
        vprint("mismatch: {}!={}?".format(pi1[worst+i],pi2[worst+j]));
        #and search for this ahead:
        recheck = np.isclose(pi1[worst+i],pi2[worst+j+1:]);#intentionally after jend.
        addedrecheck = recheck.T[0] & recheck.T[1];
        found = np.argmax(addedrecheck);
        vprint("found-j: {},{}".format(worst+found+j+1,recheck[found]));
        vprint(pi2[worst+found+j+1]);
        #if found, add particles betwen j and j+jj
        #as nans to the first *before* i, these
        #are new particles
        quit();
        if addedrecheck[found]:
            vprint("found at {}, {}".format(j+found,pi2[j+1+found]));
            vprint("adding to first which is {}+{}={}".format(j,found+1,j+found+1));
            nanrows1.append(data2[j:j+found+1]);
            nani1.append(i);
            #resume loop from next i and jj+1
            i,j = i+1, j+found+1;
        else:
            #not found, so particle lost, add it
            #to second before j
            vprint("not found, so adding to second.");
            nanrows2.append(data1[i+worst]);
            nani2.append(j+worst);
            i+=worst+1; j+=worst;
        vprint("----");
    if i<len(data1): #extra ones lost
        nanrows2.append(data1[i:])
        nani2.append(len(data2));
    elif j<len(data2): #extra new ones
        nanrows1.append(data2[j:])
        nani1.append(len(data1));

    vprint("nani2:{}".format(nani2));
    vprint("nanrows2:");
    vprint("{}".format(np.array(nanrows2)));
    return nani1,nanrows1,nani2,nanrows2
    vprint("before");
    vprint(data2.shape);

def getfills2(pi1,pi2):
    #indices along the two arrays
    j=0;
    #nan insertions, rows for data and indices
    nanrows1=[]; nani1=[];
    nanrows2=[]; nani2=[];
    def close(a,b):
        ret = np.isclose(a,b);
        return ret.T[0] and ret.T[1];
    def lessthan(a,b):
        if a[1] < b[1]: return True;
        elif a[1] > b[1]: return False;
        #iffy here, otherwise it's close? hmm...
        return a[0] < b[0];
    vprint("starting comparison");
    vprint("size of first:{}\nsize of second:{}".format(len(pi1),len(pi2)));        
    for i,pi1i in enumerate(pi1):
        if j > len(pi2): break;
        vprint("i,j={},{}".format(i,j));
        if close(pi1i, pi2[j]):
            j+=1; continue;
        else:
            jj = j+1; found = None;
            while True:
                if not lessthan(pi1i, pi2[jj]): break;
                if close(pi1i, pi2[jj]):
                    found = jj;
                    break;
            if found is None:
                #not found, so particle lost, add it
                #to second before j
                vprint("not found, so adding to second.");
                nanrows2.append(data1[i]);
                nani2.append(j);
            else:
                vprint("found at {}, {}".format(j+found,pi2[j+1+found]));
                vprint("adding to first which is {}+{}={}".format(j,found+1,j+found+1));
                nanrows1.append(data2[j:j+found+1]);
                nani1.append(i);
                j+=found+1;
        vprint("----");
    if i<len(data1): #extra ones lost
        nanrows2.append(data1[i:])
        nani2.append(len(data2));
    elif j<len(data2): #extra new ones
        nanrows1.append(data2[j:])
        nani1.append(len(data1));

    vprint("nani2:{}".format(nani2));
    vprint("nanrows2:");
    vprint("{}".format(np.array(nanrows2)));
    return nani1,nanrows1,nani2,nanrows2

def getfillshack(pi1,pi2):
    i,j=0,0;
    #nan insertions, rows for data and indices
    nanrows1=[]; nani1=[];
    nanrows2=[]; nani2=[];

    while True:
        if i > len(pi1) or j > len(pi2): break;
        iend,jend = len(pi1),len(pi2);
        if iend - i > jend - j:
            iend -= (iend - i - jend + j);
        elif iend - i < jend - j:
            jend -= (jend - j - iend+ i);
        check = np.isclose(pi1[i:iend],pi2[j:jend]);
        addedcheck = check.T[0] & check.T[1];
        vprint("i,j={},{}".format(i,j));
        worst = np.argmin(addedcheck);
        if addedcheck[worst] == True:
            i=iend;
            j=jend;
            break;
        #not found, so particle lost, add it
        #to second *after* j
        vprint("not found, so adding to second.");
        nanrows2.append(data1[i+worst]);
        nani2.append(j+worst);
        i+=worst+1; j+=worst;
        vprint("----");
    if i<len(data1): #extra ones lost
        nanrows2.append(data1[i:])
        nani2.append(len(data2));
    return nani2,nanrows2
vprint("first shape: {}".format(data1.shape));
vprint("data2 shape before: {}".format(data2.shape));
#first, simple check.

eq = lambda a,b: np.logical_and.reduce(np.isclose(a,b).T).min();
if pi1.shape != pi2.shape or not eq(pi1,pi2):
    nani2,nanrows2 = getfillshack(pi1,pi2);
    second[0]['data']=np.insert(data2, nani2, nanrows2);
else:
    vprint("equal, doing nothing");
vprint("data2 shape after: {}".format(data2.shape));
vprint("outputting...");
with open(opts['<output>'],'w') as f:
    pickle.dump(second,f,2);
pass

