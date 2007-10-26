#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# 1 - read the output
# 2 - search all floating point expressions
# 3 - replace it to have a comparable text
# 4 - compare each floating point expressions
# Date: 31/10/2006
#----------------------------------------------------------------------------

import difflib
import getopt
import re
import sys

re_float = re.compile("([-]?[0-9]+[.][0-9]+([EDed][-+]?[0-9]+)?)")

def give_text_floats(text):
    floats = list()
    new = ""
    for line in text.splitlines(1):
        if "CPU time" in line or "MB" in line or "proc" in line or "Processes" in line:
            continue
        new += line
        search = re_float.findall(line)
        for (one,two) in search:
            floats.append(one)
    #Replace all floating point by XXX
    new = re_float.sub('XXX',new).splitlines(1)
    floats = map(lambda x: float(x),floats)
    return (new,floats)

def usage():
    print "fldiff.py file1 file2"
    sys.exit(1)

if __name__ == "__main__":
    #Check arguments
    try:
        optlist, args = getopt.getopt(sys.argv[1:],[])
    except getopt.error:
        sys.stderr.write("Error in arguments\n")
        usage()
    if len(sys.argv) != 3:
        sys.stderr.write("Error in arguments\n")
        usage()
    original1 = open(sys.argv[1]).read()
    original2 = open(sys.argv[2]).read()
    (text1,floats1) = give_text_floats(original1)
    (text2,floats2) = give_text_floats(original2)
    diff = difflib.unified_diff(text1,text2) 
    print ''.join(diff)
    n = len(floats1)
    if n == len(floats2):
        discrepancy = 0.0
        for i in range(n):
            discrepancy = max(discrepancy,abs(floats1[i]-floats2[i]))
        print "Max Discrepancy: ",discrepancy
    else:
        print "The number of floating point differs"
