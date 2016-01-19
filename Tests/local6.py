# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:57:20 2016

@author: Scott
"""
import sys
import numpy as np

try:
    import AnalyzeAll as aa
except:
    print "Modifying path to include LSPreader"
    readerpath = r'..'
    sys.path.append(readerpath) # Add the LSP reader as taking precedence after all but the current directory
    import AnalyzeAll as aa

import lstools as ls
import lspreader2 as rd
    
p4root = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump'
outroot = r'C:\Users\Scott\Documents\temp\sclrtest\analysis'
#aa.analyzeAll(p4root, outroot)


p4dir = r'C:\Users\Scott\Documents\temp\sclrtest\lspdump\test0'
data = ls.readFldScl(p4dir)


chk_fs = 1.0#fs How many femtoseconds per data chunk for analysis?
data_chunked = ls.chunkData(data, chk_fs)
