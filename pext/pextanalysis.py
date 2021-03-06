# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:20:31 2016

LSP particle extraction plane analysis.

Scott's consolidated analysis of pext.p4 files (extractions of electrons from four edges.) Very customized.

@author: Scott
"""

# Matplotlib stuff
import matplotlib as mpl
#mpl.use('Agg')
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

# Other dependencies
import numpy as np
import os

import angular # Radial histogram plot, from Gregory's local file 'angular.py'
import quantities
import lspreader2 as rd # Hope this works!
import lstools as ls # local
import sftools as sf # local
import scipy.constants as sc

def pextFull(p4dir, outdir = '.', shortname = 'Sim', Utot_Jcm = None):
    """ Perform the full-blown analysis of pext*.p4(.gz) files from a simulation. 
    Inputs:
        p4dir: string, path to the folder which contains the pext*.p4
        outdir: string, path to an already-existent folder in which to save outputs (PNG plots)
        shortname: string, a fairly short identifier for the simulation which will be added into plots and filenames
        Utot_Jcm: number (optional): If total simulation energy is given (in J/cm, for a 2D sim), then relative efficiencies will be plotted and written to CSV.
    Outputs:
        Saves a bunch of well-labeled PNG plot files into the folder 'outdir'
    """
    
    fns = ls.getp4(p4dir, prefix='pext')
    if len(fns) < 1:
        msg = "No pext*.p4(.gz) files found in the directory: " + p4dir
        raise Exception(msg)
    pextarr = loadPexts(fns)
    plotme(pextarr, outdir=outdir, shortname=shortname, Utot_Jcm = Utot_Jcm) # Plots and saves
    return pextarr
    
def loadPexts(fns):
    """ Calls lspreader2 to load all the pext*.p4 from a given list of filenames, and concatenates into a single array
    Inputs:
        fns: List of filenames
    Outputs:
        pextarr: A NumPy array of the same form as the output of lspreader2.read_pext(), containing a concatenation of all the pext file contents as opposed to just one
    """
    pextarr = None # An array which will hold the contents of all the pext files
    for i in range(len(fns)):
        fn = fns[i] # Current filename
        pextarr_tmp = rd.read_pext2(fn)
        if i == 0:
            pextarr = pextarr_tmp # Initialize the output array 'outs'
        else:
            pextarr = np.concatenate((pextarr, pextarr_tmp), 0) # Append this pext*.p4 'out' to the output array 'outs'
    return pextarr

def plotme(pextarr, outdir='.', shortname='Sim', mksub=False, Utot_Jcm = None, spot3D=2.26e-4):
    """
    Make Scott's preferred set of plots (and save to file) from the output (or concatenated outputs) of lspreader "pext()" of pext.
    Inputs:
        pextarr: NumPy array, of the fancy form output by pextarr = lspreader2.read_pext('pext1.p4'). Several can be concatenated as desired.
        outdir: string, path to a folder (already existent) that can hold the outputs
        shortname: string, a fairly short identifier for the simulation which will be added into plots and filenames
        mksub: bool, if True then make a subdirectory in 'outdir', otherwise just dump the PNGs directly into 'outdir'
        Utot_Jcm: number (optional): If total simulation energy is given (in J/cm, for a 2D sim), then relative efficiencies will be plotted.
		spot3D: number, size of laser spot that would be used in a 3D simulation (cm). Gives the conversion of 2D nC/cm to 3D nC
    Outputs:
        Will save PNG figures summarizing the extraction planes into the folder specified by 'outdir'.
    
    TODO: stretch out the time of flight figure
    TODO: fix up the time of flight figure so it's not such an ugly, long string of plotting functions (copied and pasted text)
    """

    # Make a subdirectory for the PNG files, if mksub was set to True
    if mksub:
        plotdir = sf.subdir(outdir, 'Extraction planes')
    else:
        plotdir = outdir
        
    ## GET PHI FROM LSPREADER OUTPUTS OF PEXT (Concatenated into one, of course.)
    massE = 0.511 # electron rest mass, in MeV
    d = {}
    d['t'] = pextarr['t']*1e6 # Time, in fs
    u_norm = np.sqrt(pextarr['ux']**2+pextarr['uy']**2+pextarr['uz']**2) # 'uz', etc. is actually momentum / (speed of light * rest mass); it's what they call gamma-beta units; that is, ux = px/c/mass_0 = gamma * beta_x
    #p_norm = massE * u_norm # The momenta of the particles, 
    d['KE'] =(np.sqrt(u_norm**2+1)-1)*massE # Electron kinetic energy, in MeV (kinetic energy of each constituent real electron, that is; not of the whole macroparticle)
    d['q'] = -pextarr['q']*1e3 # Amount of negative charge, nanoColoumbs/cm
    d['phi'] = np.arctan2(pextarr['uz'], pextarr['ux'])
    
    m_macro = -((pextarr['q']*1e-6)/sc.e)*sc.m_e # Mass of this macroparticle, in kg/cm. Mass macroparticle = (Coulombs macroparticle / coulombs electron) * mass electron
    d['KE_macro_Jcm'] = (m_macro * sc.c**2)*(np.sqrt(u_norm**2 + 1) - 1) # Kinetic energy of this macroparticle, in J/cm
    
    if Utot_Jcm: # If the user input a value for total energy
        ## CSV: Write some basic efficiency info to file
        with open(os.path.join(outdir, shortname + ' - Electron Efficiency.csv'), 'w') as f:
            Utot_mJ3D = Utot_Jcm * spot3D * 1e3
            f.write('Total energy in sim: ' + str(Utot_Jcm) + " J/cm. (3D equivalent with spot size " + str(spot3D*1e4) + " microns is " + str(Utot_mJ3D) + " mJ.)\n")
            f.write("MINIM_ELECTRON_MEV, PERCENT_EFFICIENCY_PLUSMINUS_40DEG, PERCENT_EFFIC_PLUSMINUS_6point3_DEG, pC_3D_PLUSMINUS_40DEG, pC_3D_PLUSMINUS_6point3_DEG\n")        
            ecuts = [0.12, 0.3, 0.5, 1.0, 1.5, 3] # Cutoffs for energy efficiency, in MeV
            for ecut in ecuts:
                cdtA = np.logical_and(d['KE'] > ecut, np.abs(d['phi']) > np.deg2rad(180 - 40)) # Condition A: Wide angle
                cdtB = np.logical_and(d['KE'] > ecut, np.abs(d['phi']) > np.deg2rad(180 - 6.3)) # Condition B: Narrow angle
                Ue_JcmA = np.sum(d['KE_macro_Jcm'][cdtA]) # Energy of the electrons meeting Condition A, in J/cm
                Ue_JcmB = np.sum(d['KE_macro_Jcm'][cdtB]) # Energy of the electrons meeting Condition B, in J/cm
                efficA = Ue_JcmA/Utot_Jcm # Efficiency into condition A
                efficB = Ue_JcmB/Utot_Jcm # Efficiency into condition B
                Q_pC3DA = np.sum(d['q'][cdtA])*spot3D*1e3 # Number of electrons meeting Condition A, in pC (converted from nC/cm using spot size (cm) and factor of 1e3)
                Q_pC3DB = np.sum(d['q'][cdtA])*spot3D*1e3 # Number of electrons meeting Condition B, in pC (converted from nC/cm using spot size (cm) and factor of 1e3)
                f.write(str(ecut) + ", " + str(np.round(efficA*100,2)) + ", " + str(np.round(efficB*100,2)) + ", " + str(np.round(Q_pC3DA,2)) + ", " + str(np.round(Q_pC3DB,2)) + "\n")
            cdt120 = np.logical_and(d['KE'] > 0.120, np.abs(d['phi']) > np.deg2rad(180 - 40))
            effic120 = np.sum(d['KE_macro_Jcm'][cdt120])/Utot_Jcm
            efficstr = r'$\stackrel{>120 keV}{\pm 40^{\circ}}$ Efficiency: ' + str(np.round(effic120 * 100, 2)) + "%"
    else:
        efficstr = '' # Implying that efficiency was not calculated.
    
    ## RADIAL PLOT WITHOUT CUTOFFS
    fig = plt.figure(3)
    plt.clf()
    title=r'Radial electron spectrum'
    angular.angular2(d['q'], d['phi'], d['KE'], fig = fig, cmap='viridis', minQ=0, maxQ=1.5, title=title)
    plt.figtext(0.99, 0.02, shortname, horizontalalignment='right')
    plt.figtext(0.01, 0.02, efficstr, horizontalalignment='left')
    fig.savefig(os.path.join(plotdir, shortname + ' - Radial.png'))
   
    ## RADIAL PLOT WITH CUTOFF
    fig = plt.figure(2)
    condit2 =np.logical_and(d['KE'] > .300, np.abs(d['phi']) > np.deg2rad(180 - 40))
    plt.clf()
    #title = shortname + r": 300/40 charge: " + str(round(np.sum(d['q'][condit2])))
    title=r'Radial electron spectrum, $\pm$40$^{\circ}$'
    angular.angular2(d['q'][condit2], d['phi'][condit2], d['KE'][condit2], fig = fig, cmap='viridis', minQ=0, maxQ=1.5, title=title)
    plt.figtext(0.99, 0.02, shortname, horizontalalignment='right')
    plt.figtext(0.01, 0.02, efficstr, horizontalalignment='left')
    fig.savefig(os.path.join(plotdir, shortname + ' - Radial with cutoffs.png'))

    ## Z-ELECTRON PLOT WITH CUTOFF
    fig = plt.figure(10)
    plt.clf()
    zmax = np.max(pextarr['z']*1e4) # The time of the end of the sim, presumably. Draw a vertical line at this time.
    zmin = np.min(pextarr['z']*1e4) # The time of the end of the sim, presumably. Draw a vertical line at this time.
    zrange = (zmin,zmax) # The plotting range
    ax = plt.subplot(111)
    cdt_z = np.logical_and(np.abs(d['phi']) > np.deg2rad(180 - 40), d['KE'] > 2.00)
    [histvals, histedges] = np.histogram(pextarr['z'][cdt_z]*1e4, range=zrange, bins=800, weights=d['q'][cdt_z])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    ax.set_ylabel('Collected charge, 2 MeV min. (a.u.)')
    ax.set_xlabel('Z ($\mu m$)')
    ax.set_xlim(zrange[0], zrange[1])
    plt.title('Charge vs. Z')
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Charge vs Z.png'))


    ## Z-ORIGINATION PLOT
    try:
        pextarr['xi'] # Fail fast if this pext file doesn't have origins
        fig = plt.figure(10)
        plt.clf()
        zmax = np.max(pextarr['z']*1e4) # The time of the end of the sim, presumably. Draw a vertical line at this time.
        zmin = np.min(pextarr['z']*1e4) # The time of the end of the sim, presumably. Draw a vertical line at this time.
        zrange = (zmin,zmax) # The plotting range
        ax = plt.subplot(111)
        cdt_z = np.logical_and(np.abs(d['phi']) > np.deg2rad(180 - 40), d['KE'] > 0.12)
        [histvals, histedges] = np.histogram(pextarr['zi'][cdt_z]*1e4, range=zrange, bins=800, weights=d['q'][cdt_z])
        picYt = histvals
        picXt = histedges[1:] - (histedges[1] - histedges[0])/2
        plt.plot(picXt, picYt)
        ax.set_ylabel('Origination of charge, $>$120 keV (a.u.)')
        ax.set_xlabel('Z ($\mu m$)')
        ax.set_xlim(zrange[0], zrange[1])
        plt.title('Origination, Z slice')
        plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
        fig.savefig(os.path.join(plotdir, shortname + ' - Origination Z.png'))

    ## X-ORIGINATION PLOT
        fig = plt.figure(10)
        plt.clf()
        zmax = np.max(pextarr['x']*1e4) # The time of the end of the sim, presumably. Draw a vertical line at this time.
        zmin = np.min(pextarr['x']*1e4) # The time of the end of the sim, presumably. Draw a vertical line at this time.
        zrange = (zmin,zmax) # The plotting range
        ax = plt.subplot(111)
        cdt_z = np.logical_and(np.abs(d['phi']) > np.deg2rad(180 - 40), d['KE'] > 0.12)
        [histvals, histedges] = np.histogram(pextarr['xi'][cdt_z]*1e4, range=zrange, bins=800, weights=d['q'][cdt_z])
        picYt = histvals
        picXt = histedges[1:] - (histedges[1] - histedges[0])/2
        plt.plot(picXt, picYt)
        ax.set_ylabel('Origination of charge, $>$120 keV (a.u.)')
        ax.set_xlabel('X ($\mu m$)')
        ax.set_xlim(zrange[0], zrange[1])
        plt.title('Origination, X slice')
        plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
        fig.savefig(os.path.join(plotdir, shortname + ' - Origination X.png'))
    

    ## X/Z (2D) ORIGINATION COLORMAP
        pextarr['xi'] # Fail fast if this pext file doesn't have origins
        xdims = np.array([-8, -2]) # In microns
        zdims = np.array([-3, 3]) # In microns
        cdt_z = np.logical_and(np.abs(d['phi']) > np.deg2rad(180 - 40), d['KE'] > 0.12)
        H_now, xedges, yedges = np.histogram2d(pextarr['xi'][cdt_z]*1e4, pextarr['zi'][cdt_z]*1e4, weights=d['q'][cdt_z], range = [xdims, zdims], bins=(100,100))
    
        cmax = np.max(H_now[:])
        fig = plt.figure(11)
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_title('Origination of charge $>$120 keV (a.u.)')
        X, Y = np.meshgrid(xedges, yedges)
        im = ax.pcolorfast(X, Y, H_now.swapaxes(0,1), cmap='viridis', vmin=0, vmax=cmax)
        plt.axis('equal')
        plt.xlim(xdims)
        plt.ylim(zdims)
        plt.xlabel("X ($\mu$m)")
        plt.ylabel("Z ($\mu$m)")
        cbar = fig.colorbar(im, label="$>$120 keV Electrons/volume (a.u.)")
        #ax.text(np.min(X[:]) + 1, np.max(Y[:]) - 3, tstring, fontsize=24, color='white')
        #plt.figtext(0.01, 0.02, Qstring, horizontalalignment='left')
        plt.figtext(1 - 0.01, 0.02, shortname, horizontalalignment='right')
        fig.savefig(os.path.join(plotdir, shortname + ' - Origination 2D (120 keV).png'))
        
        # Make a 2 MeV version
        cdt_z = np.logical_and(np.abs(d['phi']) > np.deg2rad(180 - 40), d['KE'] > 1.5)
        H_now, xedges, yedges = np.histogram2d(pextarr['xi'][cdt_z]*1e4, pextarr['zi'][cdt_z]*1e4, weights=d['q'][cdt_z], range = [xdims, zdims], bins=(150,150))
    
        cmax = np.max(H_now[:])
        fig = plt.figure(12)
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_title('Origination of charge $>$1.5 MeV (a.u.)')
        X, Y = np.meshgrid(xedges, yedges)
        im = ax.pcolorfast(X, Y, H_now.swapaxes(0,1), cmap='viridis', vmin=0, vmax=cmax)
        plt.axis('equal')
        plt.xlim(xdims)
        plt.ylim(zdims)
        plt.xlabel("X ($\mu$m)")
        plt.ylabel("Z ($\mu$m)")
        cbar = fig.colorbar(im, label="$>$1.5 MeV Electrons/volume (a.u.)")
        #ax.text(np.min(X[:]) + 1, np.max(Y[:]) - 3, tstring, fontsize=24, color='white')
        #plt.figtext(0.01, 0.02, Qstring, horizontalalignment='left')
        plt.figtext(1 - 0.01, 0.02, shortname, horizontalalignment='right')
        fig.savefig(os.path.join(plotdir, shortname + ' - Origination 2D (1_5 MeV).png'))
        
    except:
        print "ORIGINATION PLOTS FAILED. Particle extraction file does not contain initial conditions (SAK_EXTRACT_IPP) (or, some other error in these plots.) Moving on to other plots."
        pass
    
    ## BACKWARD ELECTRON SPECTRUM LINE PLOT
    fig = plt.figure(4)
    plt.clf()
    cdt_inner = np.abs(d['phi']) > np.deg2rad(180 - 6.3)
    cdt_outer = np.abs(d['phi']) > np.deg2rad(180 - 40)
    nbins = 60
    Erange = np.array([0, 3.0])
    dE = (np.max(Erange) - np.min(Erange))/nbins # Width, in MeV, of each histogram bin
    [histvals, histedges] = np.histogram(d['KE'][cdt_outer], range=Erange, bins=nbins, weights=d['q'][cdt_outer])
    picY1 = histvals / dE # convert from charge/bin to nC/MeV/cm
    picX1 = histedges[1:] - (histedges[1] - histedges[0])/2
    [histvals, histedges] = np.histogram(d['KE'][cdt_inner], range=Erange, bins=nbins, weights=d['q'][cdt_inner])
    picY2 = histvals / dE # convert from charge/bin to nC/MeV/cm
    picX2 = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picX1, picY1, picX2, picY2)
    plt.xlabel("Kinetic energy, MeV")
    plt.ylabel("Charge density, 2D (nC/MeV/cm)")
    plt.ylim(0, 625)
    plt.legend([r'Backward $\pm$40$^{\circ}$', r'Backward $\pm$6.3$^{\circ}$'])
    plt.title(r'Electron spectrum', fontsize=20)
    plt.tight_layout()
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Electron spectrum.png'))

    
    ## TIME OF FLIGHT PLOT
    fig = plt.figure(5)
    tmax = np.max(d['t']) # The time of the end of the sim, presumably. Draw a vertical line at this time.
    trange = (100,400) # The plotting range
    plt.clf()
    
    ax = plt.subplot(411)
    cdt_time = np.logical_and(cdt_outer, d['KE'] > 3.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    plt.title('KE > 3 MeV')
    
    
    ax = plt.subplot(412)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 3.0, d['KE'] > 1.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    plt.title('1 MeV < KE < 3 MeV')
    
    ax=  plt.subplot(413)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 1.0, d['KE'] > .5)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    plt.title('500 keV < KE < 1 MeV')
    
    ax = plt.subplot(414)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < .5, d['KE'] > .1)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=800, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    plt.xlabel("Time (fs)")
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    plt.title('100 keV < KE < 500 keV')
    
    plt.suptitle('Electron time of flight, $\pm$40$^{\circ}$', fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Electron time of flight.png'))


    ## TIME OF FLIGHT PLOT, Zoomed in
    fig = plt.figure(6)
    tmax = np.max(d['t']) # The time of the end of the sim, presumably. Draw a vertical line at this time.
    tmid = np.sum(d['t']*d['q'])/np.sum(d['q']) # Centroid of charge deposition, in time
    trange = (tmid - 50,tmid + 50) # The plotting range
    plt.clf()
    
    ax = plt.subplot(411)
    cdt_time = np.logical_and(cdt_outer, d['KE'] > 3.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('KE > 3 MeV')
    
    
    ax = plt.subplot(412)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 3.0, d['KE'] > 1.0)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('1 MeV < KE < 3 MeV')
    
    ax=  plt.subplot(413)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < 1.0, d['KE'] > .5)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('500 keV < KE < 1 MeV')
    
    ax = plt.subplot(414)
    cdt_time = np.logical_and(cdt_outer, d['KE'] < .5, d['KE'] > .1)
    [histvals, histedges] = np.histogram(d['t'][cdt_time], range=trange, bins=300, weights=d['q'][cdt_time])
    picYt = histvals
    picXt = histedges[1:] - (histedges[1] - histedges[0])/2
    plt.plot(picXt, picYt)
    plt.axvline(tmax, color='black')
    plt.xlabel("Time (fs)")
    ax.set_xlim(trange[0], trange[1])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    ax.set_axis_bgcolor((1.0, 0.9, 0.9))
    plt.title('100 keV < KE < 500 keV')
    
    plt.suptitle('Electron time of flight (ZOOMED), $\pm$40$^{\circ}$', fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.figtext(0.99, 0.01, shortname, horizontalalignment='right')
    fig.savefig(os.path.join(plotdir, shortname + ' - Zoomed TOF.png'))

    plt.close('all') #Close all the figures we just opened to avoid hogging memory