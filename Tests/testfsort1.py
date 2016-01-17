# -*- coding: utf-8 -*-
"""
Created on Thu Jan 07 11:30:19 2016

@author: Scott
"""
# Documentation for reg exp of strings at: https://docs.python.org/2/library/re.html

# A regular expression (or RE) specifies a set of strings that matches it; the functions in this module let you check if a particular string matches a given regular expression


import re

name = 'fld5020.p4.gz'

m = re.search('c.*', 'abcdef')
if m:
    print m.group(0)
else:
    print 'No match.'
    

pattern = r'^(fld)(\d*?)(\.p4)(\.gz){0,1}$'

m = re.search(pattern, name)
# Matches 'fldXXX.p4' and 'fldXXX.p4.gz', where XXX is a number of any length

if m:
    print m.group(0)
    # Note: re.split(r'^(fld)(\d*?)(\.p4)(\.gz){0,1}$', 'fld5020.p4.gz) returns something like ['', 'fld', '5020', '.p4', '.gz', '']
    fnum = int(re.split(pattern, name)[2])
    print fnum
else:
    print 'No match.'