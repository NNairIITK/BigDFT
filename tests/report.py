#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Final report
# Date: 29/10/2007
#----------------------------------------------------------------------------

import glob
import re
import sys

#Max value
max_discrepancy = 1.e-11
#Regular expression
re_discrepancy = re.compile("Max Discrepancy:[ ]+(.*)")

print "Final report:"
for file in glob.glob("*/fldiff.report"):
    dir = file.split("/")[0]
    discrepancy = float(re_discrepancy.findall(open(file).read())[0])
    if discrepancy <= max_discrepancy:
        state = "succeeded"
    else:
        state = "failed"
    print "%-9s %s (%7.2e)" % (dir,state,discrepancy)
