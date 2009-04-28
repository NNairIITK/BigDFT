#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# 1 - read the output
# 2 - search all floating point expressions
# 3 - replace it to have a comparable text
# 4 - compare each floating point expressions
# Date: 22/04/2009
#----------------------------------------------------------------------------

import difflib
import getopt
import os
import re
import sys

#Check the version of python
version = map(int,sys.version_info[0:3])
if  version < [2,3,0]:
    sys.stderr.write("Detected version %d.%d.%d\n" % tuple(version))
    sys.stderr.write("Minimal required version is python 2.3.0: Use the command diff\n")
    os.system("diff %s %s" % (file1,file2))
    sys.exit(1)

#Match the version number ex. 1.1.9
re_version = re.compile("[(]ver[ ]+[0-9.]+[)]",re.IGNORECASE)
#Match a floating number
re_float = re.compile("([- ]?[0-9]+[.][0-9]+([EDed][-+]?[0-9]+)?)")

#Maximum discrepancy between float results (default)
max_discrepancy = 1.1e-10

def usage():
    print "fldiff.py [--bigdft] [--discrepancy=d] [--help] file1 file2"
    print "  --bigdft to compare 'bigdft' output files"
    print "  --discrepancy=%7.1e Maximal discrepancy between results" % max_discrepancy
    print "  --help   display this message"
    sys.exit(1)

#Maximum discrepancy between float results
max_discrepancy = 1.1e-10

#Check arguments
try:
    optlist, args = getopt.getopt(sys.argv[1:],"bd:h",["bigdft","discrepancy=","help"])
except getopt.error:
    sys.stderr.write("Error in arguments\n")
    usage()
bigdft = False
for opt,arg in optlist:
    if opt == "-b" or opt == "--bigdft":
        bigdft = True
    elif opt == "-d" or opt == "--discrepancy":
        max_discrepancy=float(arg)
    elif opt == "-h" or opt == "--help":
        usage()
if len(args) != 2:
    sys.stderr.write("Error in arguments\n")
    usage()

#Display the maximum discrepancy
print max_discrepancy

#Arguments
file1 = args[0]
file2 = args[1]

if bigdft:
    #Test if the line should not be compared (bigdft output)
    def line_junk(line):
        "True if the line must not be compared"
        return re_version.search(line) \
            or "CPU time" in line \
            or "Load" in line \
            or "memory" in line \
            or "MB" in line \
            or "proc" in line \
            or "Processes" in line \
            or "allocation" in line \
            or "~W" in line \
            or "for the array" in line
else:
    def line_junk(line):
        "Always False"
        return False

#Check the last line
end_line = "MEMORY CONSUMPTION REPORT"

#Read the first file
try:
    original1 = open(file1).read().replace('\r','').splitlines(1)
except IOError:
    sys.stderr.write("The file '%s' does not exist!\n" % file1)
    sys.exit(1)
#Read the second file
try:
    original2 = open(file2).read().replace('\r','').splitlines(1)
except IOError:
    sys.stderr.write("The file '%s' does not exist!\n" % file2)
    sys.exit(1)

maximum = 0.0
context_discrepancy = ""
context_lines = ""

#First we compare the first two lines in the case of an added prefix 
#(as in platine computer of CCRT)
#We detect a pattern
if bigdft:
    pattern = '                             BBBB         i       ggggg    '
else:
    pattern = ''

try:
    p1 = original1[0].index(pattern)
    p2 = original2[0].index(pattern)
except ValueError:
    #First line not found ??
    p1 = -1
    p2 = -1
except IndexError:
    sys.stdout.write("One file is blank!\n")
    sys.exit(1)

if p1 >= 0 and p2 >= 0 and p1 != p2:
    #we try something
    prefix1 = original1[0].replace(pattern,'')[:-1]
    prefix2 = original2[0].replace(pattern,'')[:-1]
    if prefix1 != '':
        #We remove the prefix
        original1 = map(lambda x: x.replace(prefix1,''), original1)
    if prefix2 != '':
        #We remove the prefix
        original2 = map(lambda x: x.replace(prefix2,''), original2)

end_left = False
for line in original1:
    end_left = end_line in line
    if end_left:
        break
end_right = False
for line in original1:
    end_right = end_line in line
    if end_right:
        break

if bigdft:
    #Do not compare if a file is not properly finished
    if not end_left:
        print "WARNING: The file '%s' is not properly finished!" % file1
    if not end_right:
        print "WARNING: The file '%s' is not properly finished!" % file2
    if not (end_left and end_right): 
        print "Max Discrepancy: NaN"
        sys.exit(1)

#Compare both files
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
    print_context = False
    left = list()
    line = compare.next()
    while line[0] == "-":
        left.append(line)
        try:
            line = compare.next()
        except StopIteration:
            #We have reached the end of file
            EOF = True
            break
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
        if line_junk(line1) or line_junk(line2):
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
        if new1 != new2 and i1 == 0:
            #For the first difference, we display the context
            print context,
        print_context = True
        if new1 != new2:
            print line1,
            print line2,
        n = len(floats1)
        if n == len(floats2):
            diff_discrepancy = False
            for i in range(n):
                tt = abs(floats1[i]-floats2[i])
                if maximum < tt:
                    context_discrepancy = " (line %s)" % context.split(",")[0][4:]
                    context_lines = "\n"+context_discrepancy[1:]+"\n"+line1+line2
                    maximum = max(maximum,tt)
                if tt > max_discrepancy:
                    diff_discrepancy = True
            if diff_discrepancy and new1 == new2:
                if not print_context:
                    print context,
                    print_context = True
                print line1,
                print line2,
        else:
            print "%s the number of floating point differs" % context[:-1]
    #Add lines if necessary
    while i1 < n1-1:
        i1 += 1
        if n1 > 0 and not line_junk(left[i1]):
            if not print_context:
                print context,
            print_context = True
            print left[i1],
    while i2 < n2-1:
        i2 += 1
        if n2 > 0 and not line_junk(right[i2]):
            if not print_context:
                print context,
            print_context = True
            print right[i2],

print context_lines,

if maximum > max_discrepancy:
    start = "\033[0;31m"
    message = "failed    < "
    end = "\033[m"
else:
    start = "\033[0;32m"
    message = "succeeded < "
    end = "\033[m"

print "%sMax Discrepancy %s: %s (%s%s)%s" % (start,context_discrepancy,maximum,message,max_discrepancy,end)
sys.exit(0)

