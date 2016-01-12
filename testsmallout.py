from h5stitch2D import *
from freqanalysis import *
import lspreader2 as rd
import numpy as np
import os


# Matplotlib stuff
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from distutils.version import LooseVersion
if LooseVersion(mpl.__version__) < LooseVersion('1.5.0'):    
    # http://stackoverflow.com/questions/11887762/how-to-compare-version-style-strings and 
    print "Matplotlib", mpl.__version__, "might not have colormap 'viridis'. Importing from local colormaps.py."
    import colormaps as cmaps
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    plt.register_cmap(name='inferno', cmap=cmaps.inferno)
    plt.register_cmap(name='magma', cmap=cmaps.magma)
    plt.register_cmap(name='plasma', cmap=cmaps.plasma)




#p4dir = r'/media/sf_temp/smallout'
p4dir = r'/media/sf_temp/lastest/run2_lastest longer/gzipdir/test'
pltdir = r'/media/sf_temp/smallout/plots'

fns = getfnsp4(p4dir)
freqFull(p4dir, outdir = pltdir)
sys.exit('All done.')
