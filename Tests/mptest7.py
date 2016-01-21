import sys
try:
    import lspreader2 as rd
except:
    print "Modifying path to include LSPreader"
    readerpath = r'/home/feister.7/lsp/lspreader'
    sys.path.append(readerpath) # Add the LSP reader as taking precedence after all but the current directory
    import lspreader2 as rd
import lstools as ls
import numpy as np
import multiprocessing

if __name__=='__main__':
	p4dir = r'/data/feister.7/lspdump/a50f-14xL_mres_so-2016-01-17_2006'
	fns = ls.getp4(p4dir)
	print len(fns)

	p = multiprocessing.Pool(10)
	ls.fields2D(fns[0:1400], pool=p)
	#ls.fields2D(fns[0:700])
