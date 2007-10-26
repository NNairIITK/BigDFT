#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# 1 - read the output
# 2 - search all floating point expressions
# 3 - replace it to have a comparable text
# 4 - compare each floating point expressions
# Date: 26/10/2006
#----------------------------------------------------------------------------

import difflib
import getopt
import re
import sys

re_float = re.compile("([-]?[0-9]+[.][0-9]+([EDed][-+]?[0-9]+)?)")

def give_text_floats(text):
    return (new,floats)

def usage():
    print "fldiff.py file1 file2"
    sys.exit(1)

#Check arguments
try:
    optlist, args = getopt.getopt(sys.argv[1:],[])
except getopt.error:
    sys.stderr.write("Error in arguments\n")
    usage()
if len(sys.argv) != 3:
    sys.stderr.write("Error in arguments\n")
    usage()
#Read The first file
original1 = open(sys.argv[1]).read().splitlines(1)
#Read the second file
original2 = open(sys.argv[2]).read().splitlines(1)
n1 = len(original1)
n2 = len(original2)
i1 = 0
i2 = 0
discrepancy = 0.0
while i1 < n1:
    line1 = original1[i1]
    line2 = original2[i2]
    if "CPU time" in line1 or "MB" in line1 or "proc" in line1 or "Processes" in line1:
        pass
    else:
        floats1 = list()
        for (one,two) in re_float.findall(line1):
            floats1.append(float(one))
        floats2 = list()
        for (one,two) in re_float.findall(line2):
            floats2.append(float(one))
        #Replace all floating point by XXX
        new1 = re_float.sub('XXX',line1)
        new2 = re_float.sub('XXX',line2)
        if new1 != new2:
            print "l",i1
            print "< ", line1,
            print "> ", line2,
        n = len(floats1)
        if n == len(floats2):
            for i in range(n):
                discrepancy = max(discrepancy,abs(floats1[i]-floats2[i]))
        else:
            print "The number of floating point differs"
    i1 += 1
    i2 += 1

print "Max Discrepancy: ",discrepancy
