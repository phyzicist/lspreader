# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 13:17:24 2016

@author: Scott
"""
import pandas as pd
import numpy as np

fn = r'C:\Users\Scott\Box Sync\Analysis3\a28f-10_mres_so\a28f-10_mres_so - Electron Efficiency.csv'
dat1 = pd.read_csv(fn, skiprows=1, index_col=0)

fn = r'C:\Users\Scott\Box Sync\Analysis3\a28f-10_mres_so\a28f-10_mres_so - Reflectivity.csv'
dat2 = pd.read_csv(fn, skiprows=0, index_col=0)

print dat1, dat2