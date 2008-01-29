#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Final report
# Date: 29/01/2008
#----------------------------------------------------------------------------

import glob
import re
import sys

#Max value
max_discrepancy = 5.0e-10
#Regular expression
re_discrepancy = re.compile("Max Discrepancy[^:]*:[ ]+(.*)")

Exit = 0
print "Final report (max discrepancy=%7.1e):" % max_discrepancy
for file in glob.glob("*/fldiff.report"):
    dir = file.split("/")[0]
    discrepancy = re_discrepancy.findall(open(file).read())
    if discrepancy:
        discrepancy = float(discrepancy[0])
        if discrepancy <= max_discrepancy:
            state = "succeeded"
        else:
            state = "failed"
            Exit = 1
        print "%-9s %-9s (%7.1e)" % (dir,state,discrepancy)
#Error code
sys.exit(Exit)
