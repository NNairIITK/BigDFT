#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Final report (read fldiff.report from fldiff.py)
# Date: 25/12/2008
#----------------------------------------------------------------------------

import glob
import re
import sys

#Regular expression
re_discrepancy = re.compile("Max Discrepancy[^:]*:[ ]+([^ ]+)")

Exit = 0
print "Final report:"
for file in glob.glob("*/fldiff.report"):
    dir = file.split("/")[0]
    #Max value
    max_discrepancy = float(open(file).readline())
    discrepancy = re_discrepancy.findall(open(file).read())
    if discrepancy:
        discrepancy = float(discrepancy[0])
        if discrepancy <= max_discrepancy:
            start = "\033[0;32m"
            state = "succeeded < %7.1e" % max_discrepancy
            end = "\033[m"
        else:
            start = "\033[0;31m"
            state = "failed    > %7.1e" % max_discrepancy
            end = "\033[m"
            Exit = 1
        print "%s%-13s %-9s (%7.1e)%s" % (start,dir,state,discrepancy,end)
#Error code
sys.exit(Exit)
