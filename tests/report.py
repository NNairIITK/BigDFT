#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Build the final report (read *.report from fldiff.py)
# Date: 11/09/2012
#----------------------------------------------------------------------------

import fnmatch
import os
import re
import sys
import yaml

#Regular expressions
re_discrepancy = re.compile("Max [dD]iscrepancy[^:]*:[ ]+([^ ]+)[ ]+\(([^ ]+)")
re_time = re.compile("-- time[ ]+([0-9.]+)")

def callback(pattern,dirname,names):
    "Return the files given by the pattern"
    for name in names:
        if fnmatch.fnmatch(name,pattern):
            files.append(os.path.join(dirname,name))

def yaml_callback(pattern,dirname,names):
    "Return the files given by the pattern"
    for name in names:
        if fnmatch.fnmatch(name,pattern):
            yaml_files.append(os.path.join(dirname,name))


#List of files
files = []
os.path.walk(".",callback,"*.report")
#Sort files
files.sort()

#List of files (yaml case)
yaml_files = []
os.path.walk(".",yaml_callback,"*.report.yaml")
#Sort files
yaml_files.sort()

#Check if the output is a tty to print in colour
if sys.stdout.isatty():
    start_fail = "\033[0;31m"
    start_success = "\033[0;32m"
    start_pass = "\033[0;33m"
    end = "\033[m"
else:
    start_fail = ""
    start_success = ""
    start_pass = ""
    end = ""

#Error code
Exit = 0

#Total time for the tests
totime=0

print "Final report for writings in stdout ('passed' means all significant floats are correct):"
for file in files:
    dir = os.path.normpath(os.path.dirname(file))
    fic = "(%s)" % os.path.basename(file)
    #Max value
    try:
        max_discrepancy = float(open(file).readline())
        line = open(file).read()
        discrepancy = re_discrepancy.findall(line)
    except:
        discrepancy = False
        try:
            should_compare=open(file).readline()
            empty_file=("Nothing to compare" in should_compare)
        except:
            empty_file=False
    if discrepancy:
        #If nan gives nan (not a number and all comparisons are false)
        diff = float(discrepancy[0][0])
        if diff < max_discrepancy:
            #Four cases: 
            if discrepancy[0][1] == "failed-memory":
                #The test is OK but the memory remaining is not 0
                start = start_fail
                state = "remaining memory != 0B failed"
                Exit = 1
            elif discrepancy[0][1] == "passed":
                #passed: significant numbers (more than 5 digits) are < max_discrepancy
                start = start_pass
                state = "%7.1e < (%7.1e)    passed" % (diff,max_discrepancy)
            elif discrepancy[0][1] == "passed":
                #passed: significant numbers (more than 5 digits) are < max_discrepancy
                start = start_pass
                state = "%7.1e < (%7.1e)    passed" % (diff,max_discrepancy)
            else:
                #All numbers even with only 5 digits or less
                start = start_success
                state = "%7.1e < (%7.1e) succeeded" % (diff,max_discrepancy)
        else:
            start = start_fail
            state = "%7.1e > (%7.1e)    failed" % (diff,max_discrepancy)
            Exit = 1
        #Test if time is present
        time = re_time.findall(line)
        if time:
            totime += float(time[0])
            time = "%10ss" % time[0]
        else:
            time = ""
        print "%s%-27s %-37s %s%s%s" % (start,dir,fic,state,time,end)
    else:
        if not empty_file:
            start = start_fail
            state = "can not parse file.    failed"
            print "%s%-27s %-37s %s%s" % (start,dir,fic,state,end)

print "Final report for yaml outputs: if succeeded %53s" % "max diff (significant epsilon)"
for file in yaml_files:
    dirc = os.path.normpath(os.path.dirname(file))
    fic = "(%s)" % os.path.basename(file)
    dirfic = ("%-27s %-38s" % (dirc,fic)).strip()
    documents=[a for a in yaml.load_all(open(file, "r"), Loader = yaml.CLoader)]
    #find whether all the tests have passed (look at last part)
    try:
        discrepancy=documents[-1]["Test succeeded"]
        #test failes
        if not discrepancy:
            Exit = 1
            start = start_fail
            state = "Failed:    %7.1e > %7.1e (%s)" % \
                    (documents[-1]["Maximum discrepancy"], \
                     documents[-1]["Maximum tolerance applied"], \
                     documents[-1]["Failure reason"])
        else:
            start = start_success
            state = "Succeeded: %7.1e (%7.1e) " % \
                    (documents[-1]["Maximum discrepancy"], \
                     documents[-1]["Maximum tolerance applied"])
        #Test if time is present
        time = documents[-1]["Seconds needed for the test"]
        totime += time
        print "%s%-66s %s%8.2fs%s" % (start,dirfic,state,time,end)
    except:
        start = start_fail
        state = "Failed: Can not parse file!"
        print "%s%-66s %s%s" % (start,dirfic,state,end)


#Hours, minutes and seconds
totimeh = int(totime/3600)
totimem = int(totime-totimeh*3600)/60
totimes = totime-totimem*60-totimeh*3600
p_time  = "%sh %sm %ss" % (totimeh,totimem,totimes)
print 105*"-"
print 63*" "+"Time Needed for timed tests:%14s%s" % (p_time,end)
if Exit==0:
    print "Test set succeeded!"
else:
    print "Test set failed, check the above report!"
#Error code
sys.exit(Exit)
