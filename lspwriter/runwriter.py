# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 12:14:32 2016

Customized .lsp, .pbs, etc. run writer. Change a subset of LSP parameters in a simple way that makes sense to me.
Change the settings in the 'if __name__' section, then call this script as main.

@author: Scott

Possible issues:
* Is it okay that it outputs "8e5" instead of "8.e5", for example, or will all hell break loose?

"""

import re
import numpy as np
import scipy.constants as sc
from plaswriter import simpleColumn, columnAnalyze
import os
import shutil

def subdir(folder, name):
    """ Make a subdirectory in the specified folder, if it doesn't already exist"""
    subpath = os.path.join(folder,name)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    return subpath
    
def fileSub(template, outfile, dictionary):
    """
    MODIFY THE TEXT OF TEMPLATE FILE WITH DEFINITIONS IN DICTIONARY, CHECK FOR ERRORS, AND SAVE A NEW OUTPUT FILE
    INPUTS:
        template: string, path to template file (e.g. 'template.lsp')
        outfile: string, path to a (yet uncreated) output file (e.g. 'myrun.lsp')
        dictionary: dict, defines substitutions to make in template file (E.g. if dictionary['FLD_DUMP'] = 'ON', then in template.lsp, anywhere you see <FLD_DUMP> will become 'ON')
    OUTPUTS:
        Outputs a file to path "outfile" which contains a modified version of template.lsp (with dictionary variables substituted in. (E.g. d['FLD_DUMP'] = 'ON', then <FLD_DUMP> will become 'ON'.)
    """
    # Read in the template file
    with open(template, 'r') as f:
        s = f.read()
    
    # Search the template.lsp file, replacing the "<" + dictkey + ">" keywords with their values.
    for k in dictionary.keys():
        var = "<" + k + ">"
        if s.find(var) >= 0: 
            s = s.replace(var, dictionary[k])
        else:
            msg = "Variable " + var + " not found in template file. Either delete it from your python script, or add it somewhere the template template. Aborting."
            raise NameError(msg)
    
    # Check for errors in what was done.
    pattern = r'\<.*?\>'
    m = re.search(pattern, s)
    if m:
        msg = r"Variable " + m.group(0) + r" not found in python script. Either add it to your python script, or remove it from your template file. Aborting."
        raise NameError(msg)
    
    print "Writing dictionary-substituted file."
    with open(outfile, 'w') as f:
        f.write(s)
    
if __name__ == "__main__":
    ## USER, DEFINE THE SIMULATION WITH HIGH-LEVEL VARIABLES (Units are microns, nm, etc.)
    shortname = 'a50fp10_lres_so' # NO SPACES/slashes ALLOWED! Short name for simulation
    title = 'Hotwater in 2D I = 3e18 W cm-2, 5.0um scale, focus X=+10um, lam/8' # Simulation title
    
    scale = 5.0#um # Exponential scale length of pre-plasma
    focx = 10#um # X position of best laser focus, in microns

    # Simulation temporal
    t_f = 400#fs # maximum time, in fs
    tres = 8 # Laser cycles per timestep 
    skipt = 1 # time skip interval for field/scalar dumps. First one is always dumped.
    
    # Grid spatial
    xdims = (-35, 5)#um # X limits to simulation grid space, in microns
    xres = 8 # Laser wavelengths per cell, in X direction
    skipx = 2 # X skip interval for field/scalar dumps
    zdims = (-20, 20)#um # Z limits to simulation grid space, in microns
    zres = 8 # Laser wavelengths per cell, in Z direction 
    skipz = 2 # Z skip interval for field/scalar dumps
    
    # Pre-plasma spatial/density
    pc_xdims = (-30, 0)#um # X limits to simulation particle creation space, in microns
    pc_zdims = (-15, 15)#um # Z limits to simulation particle creation space, in microns
    columndat = 'watercolumn.dat' # filename of .dat defining electron density vs. X (file will be created below)
    
    # Laser
    sinedat = 'sine700points.dat' # filename of .dat defining laser oscillations through time, which will be copied from this folder into the output folder
    fwhm = 30#fs # Gaussian temporal FWHM of laser, in fs
    Emax = 4.763e7 # Peak electric field of laser, in kV/cm
    wlen = 800#nm # Laser wavelength, in nm
    spot = 2.26#um # Laser spot size, in microns

    # P4 outputs:
    fld = True # True if flds*.p4 outputs desired
    scl = True # True if scl*.p4 outputs desired
    pmov = False # True if pmov*.p4 outputs desired (note: these are big)
    pmov_step = 0.2#fs # Pmovie timestep interval, in fs. Only applied if pmov = True.
        
    # Supercomputer
    nodes = 1 # Number of nodes
    ppn = 48 # Number of processors per node

    ## SIMPLE CALCULATIONS
    xcells = np.round((np.max(xdims) - np.min(xdims))*(xres*1.05)/(wlen * 1e-3)) # Number of cells in the X direction. Multiply resolution by 1.05 for good measure.
    zcells = np.round((np.max(zdims) - np.min(zdims))*(zres*1.05)/(wlen * 1e-3)) # Number of cells in the Z direction. Multiply resolution by 1.05 for good measure.
    print "Number of cells in X, Z:", xcells, zcells
    dt = np.round(((wlen*1e-9)/(sc.c*1e-15))/(tres*1.05),3) # time step, in fs. (Calculate time per cycle using speed of light, then divide by desired time steps per cycle. Multiply resolution by 1.05 and round off answer to three decimal places for good measure.)
    print "Time step is:", dt, "fs"
        
    ## BUILD PBS SUBSTITUTION DICT
    pb = {}
    pb['SIMNAME'] = shortname
    pb['NODES'] = str(int(nodes)) # Number of nodes allocated
    pb['PPN'] = str(int(ppn)) # Number of nodes allocated
    pb['HOURS'] = str(int(168)) # 1-week (168 hours) maximum sim length
    pb['NPROCS'] = str(int(nodes*ppn)) # Number of processors on which to run the sim
    pb['HOSTOPT'] = '' # An added flag used only when sim spans multiple nodes    
    #pb['HOSTOPT'] = '--hostfile $MYHOSTS ' if (nodes>1) else '' # An added flag used only when sim spans multiple nodes

    ## BUILD LSP SUBSTITUTION DICT
    # a python dictionary of strings defining things in LSP language
    d = {}
    d['TITLE'] = title
    d['TIMELIM_NS'] = str(t_f*1e-6)
    d['TIMESTEP_NS'] = str(dt*1e-6)
    
    # Outputs
    d['FLDDUMP'] = 'ON' if fld else 'OFF' # Enable or disable field outputs
    d['SCLDUMP'] = 'ON' if scl else 'OFF' # Enable or disable scalar outputs, notably density
    d['PMOV_NS'] = str(pmov_step*1e-6) if pmov else '1.e+9' # particle_movie_interval_ns. 1.e+9 means off.
    d['SKIP_T'] = str(int(skipt)) # Field/scalar temporal skip in time steps
    d['SKIP_X'] = str(int(skipx)) # Field/scalar spatial skip in X
    d['SKIP_Z'] = str(int(skipz)) # Field/scalar sptial skip in Z
    
    # Sim and particle creation boundaries
    d['XMIN_CM'] = str(np.min(xdims)*1e-4)
    d['XMAX_CM'] = str(np.max(xdims)*1e-4)
    d['XCELLS'] = str(int(xcells))
    d['ZMIN_CM'] = str(np.min(zdims)*1e-4)
    d['ZMAX_CM'] = str(np.max(zdims)*1e-4)
    d['ZCELLS'] = str(int(zcells))
    
    d['PC_XMIN_CM'] = str(np.min(pc_xdims)*1e-4)
    d['PC_XMAX_CM'] = str(np.max(pc_xdims)*1e-4)
    d['PC_ZMIN_CM'] = str(np.min(pc_zdims)*1e-4)
    d['PC_ZMAX_CM'] = str(np.max(pc_zdims)*1e-4)
    d['NDOMS'] = str(int(nodes*ppn)) # Set the number of domains to the number of processors running the sim
    
    # Laser focal depth
    d['FOCX'] = str(focx*1e-4)
    
    # Laser temporal function
    d['LAS_DAT'] = sinedat # data_file (filename)
    d['LAS_INDEP'] = str(2*fwhm*1e-6) # independent_variable_multiplier (2x temporal FWHM pulse)
    d['LAS_DEP'] = str(Emax) # dependent_variable_multiplier (Emax in kV/cm units, 4763e7 => 3 x 10^18 W/cm2)
    d['WLEN_CM'] = str(wlen*1e-7)
    d['SPOT_CM'] = str(spot*1e-4)
    
    # Plasma density function
    d['DENS_DAT'] = columndat
    
    
    # Allocate filenames
    outdir = r'C:\Users\Scott\Documents\LSP\LSPwrite submissions'
    lspfn = os.path.join(subdir(outdir,shortname), shortname + '.lsp')
    pbsfn = os.path.join(subdir(outdir,shortname), 'submit.pbs')
    columnfn = os.path.join(subdir(outdir,shortname), columndat)
    sinefn = os.path.join(subdir(outdir, shortname), sinedat)
    ## MODIFY THE TEXT OF 'TEMPLATE.LSP' WITH ABOVE DEFINITIONS, CHECK FOR ERRORS, AND SAVE A NEW .LSP FILE
    fileSub('template.lsp', lspfn, d)
    fileSub('template.pbs', pbsfn, pb)
    shutil.copy(sinedat, sinefn)
    ## WRITE THE WATER COLUMN FILE
    X,Y = simpleColumn(columnfn, pc_xdims, scale = scale, npoints = 3000) # Writes to file.
    columnAnalyze(columnfn, plot=False) # Print some analysis, such as critical density X values