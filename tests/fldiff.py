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

re_float = re.compile("([- ]?[0-9]+[.][0-9]+([EDed][-+]?[0-9]+)?)")

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

file1 = sys.argv[1]
file2 = sys.argv[2]
#Read the first file
original1 = open(file1).read().splitlines(1)
#Read the second file
original2 = open(file2).read().splitlines(1)

max_discrepancy = 0.0
context_discrepancy = ""

compare = difflib.unified_diff(original1,original2,n=0)

try:
    #The 2 first lines are irrelevant
    compare.next()
    compare.next()
    line = compare.next()
    EOF = False
except StopIteration:
    #Nothing to compare
    EOF = True

while not EOF:
    if line[0] != "@":
        #Trouble
        print "Trouble:",line,
        sys.exit(1)
    #A new context is detected
    context = line
    left = list()
    line = compare.next()
    while line[0] == "-":
        left.append(line)
        line = compare.next()
    right = list()
    while line[0] == "+":
        right.append(line)
        try:
            line = compare.next()
        except StopIteration:
            #We have reached the end of file
            EOF = True
            break
    #We have a new context and two lists to compare
    n1 = len(left)
    n2 = len(right)
    i1 = -1
    i2 = -1
    while i1 < n1-1 and i2 < n2-1:
        i1 += 1
        i2 += 1
        line1 = left[i1]
        line2 = right[i2]
        #We avoid some lines
        if "CPU time" in line1 or "MB" in line1 or "proc" in line1 or "Processes" in line1 \
                or "allocation" in line1 or "~W" in line1 or "for the array" in line1:
            continue
        floats1 = list()
        for (one,two) in re_float.findall(line1):
            floats1.append(float(one))
        floats2 = list()
        for (one,two) in re_float.findall(line2):
            floats2.append(float(one))
        #Replace all floating point by XXX
        new1 = re_float.sub('XXX',line1[2:])
        new2 = re_float.sub('XXX',line2[2:])
        if new1 != new2:
            print context,
            print line1,
            print line2,
        n = len(floats1)
        if n == len(floats2):
            for i in range(n):
                tt = abs(floats1[i]-floats2[i])
                if max_discrepancy < tt:
                    context_discrepancy = " (line %s)" % context.split(",")[0][4:]
                    max_discrepancy = max(max_discrepancy,tt)
        else:
            print "The number of floating point differs"
#
#if n1 > n2:
#    print "The file '%s' has more lines (%d) than the file '%s'." % (file1,n1-n2,file2) 
#elif n2 > n1:
#    print "The file '%s' has more lines (%d) than the file '%s'." % (file2,n2-n1,file1) 
#
print "Max Discrepancy%s: " % context_discrepancy,max_discrepancy
