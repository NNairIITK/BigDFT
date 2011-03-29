#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Build the final report (read *.report from fldiff.py)
# Date: 28/03/2011
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

#List of files
files = []
os.path.walk(".",callback,"*.report")

#Check if the output is a tty to print in colour
if sys.stdout.isatty():
    start_fail = "\033[0;31m"
    start_success = "\033[0;32m"
    end = "\033[m"
else:
    start_fail = ""
    start_success = ""
    end = ""

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
            start = start_success
            state = "%7.1e < (%7.1e) succeeded" % (discrepancy,max_discrepancy)
        else:
            start = start_fail
            state = "%7.1e > (%7.1e)    failed" % (discrepancy,max_discrepancy)
            Exit = 1
        print "%s%-23s %-28s %s%s" % (start,dir,fic,state,end)
#Error code
sys.exit(Exit)
