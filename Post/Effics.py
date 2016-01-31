# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 13:17:24 2016

@author: Scott
"""
import pandas as pd
import numpy as np
import seaborn as sns
import os

analydir = r'C:\Users\Scott\Box Sync\Analysis3'

def getItems(shortname): # Relies on analydir being defined above as a global var
    fn = os.path.join(analydir, shortname, shortname + ' - Electron Efficiency.csv')
    effic = pd.read_csv(fn, skiprows=1, index_col=0)
    
    fn = os.path.join(analydir, shortname, shortname + ' - Reflectivity.csv')
    reflec = pd.read_csv(fn, skiprows=0, index_col=0)
    return effic, reflec # dataframes


shortname = r'a28f-10_mres_so'
dat1, _ = getItems(shortname)

shortname = r'a28f-26_mres_so'
dat2, _ = getItems(shortname)

frames = [dat1, dat2]
result = pd.concat(frames, keys=['x', 'y'])
