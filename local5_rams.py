#!/usr/bin/python

# Call via: "python local5_rams.py '2d-nosolid_focp5-'"

from freqanalysis import freqFull, subdir
import os
import sys

#p4dir = r'/tmp/ngirmang.1-2d-nosolid-151113/' # This folder contains the p4 files
#p4dir = /data/feister.7/lspdump/2d-nosolid-
#outdir = r'/data/feister.7/lsp/runs/2d-nosolid_foc-10' # This folder already exists and will store the files in subfolders

name = sys.argv[1]

analydir = r'/data/feister.7/lsp/analysis'
dumpdir = r'/data/feister.7/lspdump'

p4dir = os.path.join(dumpdir, name)
outdir = subdir(analydir, name)

freqFull(p4dir, outdir, nbatch = 80, divsp = 1, npool = 1)
