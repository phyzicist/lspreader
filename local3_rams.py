from freqanalysis import freqFull

p4dir = r'/tmp/ngirmang.1-2d-nosolid-151113/' # This folder contains the p4 files
outdir = r'/data/feister.7/lsp/runs/greg_run2' # This folder already exists and will store the files in subfolders

freqFull(p4dir, outdir, nbatch = 70, divsp = 1, npool = 1)