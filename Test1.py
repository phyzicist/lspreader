# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 15:30:36 2015

@author: Scott
"""

# A test of the lspreader

#import h5py
import lspreader2 as rd

fn = r"C:\Users\Scott\Documents\temp\lsp laser test fields\flds1330.p4"

data, header, doms2 = rd.read_flds2(fn, flds=['E','B'])
