# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:11:14 2016

@author: Scott
"""

import shutil
import os

def sname(foc):
    if foc <= 0:
        return 'a28f' + str(int(foc)) + '_mres_so'
    else:
        return 'a28fp' + str(int(foc)) + '_mres_so'

focs = [42, 30, -2, -6, -10, -14, -18, -22, -26, -30, -34]

analydir = r'C:\Users\Scott\Box Sync\Analysis3'
slicedir = r'C:\Users\Scott\Box Sync\Reslices\a28_midres'

nslices = len(focs)
#shortnames = [None]*nslices

for i in range(nslices):
    foc = focs[i]
    shortname = sname(foc)
    shutil.copy(os.path.join(analydir, shortname, shortname + ' - BScat Spectrum.png'), os.path.join(slicedir, 'Bscat', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, shortname + ' - Radial.png'), os.path.join(slicedir, 'Radial', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, shortname + ' - Electron time of flight.png'), os.path.join(slicedir, 'TOF', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, shortname + ' - Electron spectrum.png'), os.path.join(slicedir, 'Espec', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, shortname + ' - BScat vs Z.png'), os.path.join(slicedir, 'ZBscat', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, 'EM Omega', '01604.png'), os.path.join(slicedir, 'Omega01604', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, 'EM Omega', '02004.png'), os.path.join(slicedir, 'Omega02004', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, 'EM Omega', '01002.png'), os.path.join(slicedir, 'Omega01002', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, 'Oxygen ionization', '01103.png'), os.path.join(slicedir, 'Oxyg01103', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, 'EM Static', '01604.png'), os.path.join(slicedir, 'Static01604', str(foc+100) + '.png'))
    shutil.copy(os.path.join(analydir, shortname, 'EM Half', '01002.png'), os.path.join(slicedir, 'Half01002', str(foc+100) + '.png'))
