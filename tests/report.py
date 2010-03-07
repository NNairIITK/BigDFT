#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Final report (read fldiff.report from fldiff.py)
# Date: 06/03/2010
#----------------------------------------------------------------------------

import glob
import re
import sys

#Regular expression
re_discrepancy = re.compile("Max Discrepancy[^:]*:[ ]+([^ ]+)")

Exit = 0
print "Final report:"
for file in glob.glob("*/*.report"):
    ll = file.split("/")
    dir = ll[0]
    fic = "("+ll[1]+")"
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
        print "%s%-13s %-17s %s%s" % (start,dir,fic,state,end)
#Error code
sys.exit(Exit)
