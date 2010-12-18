#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Build the final report (read *.report from fldiff.py)
# Date: 18/12/2010
#----------------------------------------------------------------------------

import fnmatch
import os
import re
import sys

#Regular expression
re_discrepancy = re.compile("Max Discrepancy[^:]*:[ ]+([^ ]+)")

def callback(pattern,dirname,names):
    "Return the files given by the pattern"
    for name in names:
        if fnmatch.fnmatch(name,pattern):
            files.append(os.path.join(dirname,name))

#Liste of files
files = []
os.path.walk(".",callback,"*.report")

Exit = 0
print "Final report:"
for file in files:
    dir = os.path.normpath(os.path.dirname(file))
    fic = "(%s)" % os.path.basename(file)
    #Max value
    max_discrepancy = float(open(file).readline())
    discrepancy = re_discrepancy.findall(open(file).read())
    if discrepancy:
        discrepancy = float(discrepancy[0])
        if discrepancy <= max_discrepancy:
            start = "\033[0;32m"
            state = "%7.1e < (%7.1e) succeeded" % (discrepancy,max_discrepancy)
            end = "\033[m"
        else:
            start = "\033[0;31m"
            state = "%7.1e > (%7.1e)    failed" % (discrepancy,max_discrepancy)
            end = "\033[m"
            Exit = 1
        print "%s%-20s %-28s %s%s" % (start,dir,fic,state,end)
#Error code
sys.exit(Exit)
