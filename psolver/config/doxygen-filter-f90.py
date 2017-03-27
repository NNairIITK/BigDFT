#!/usr/bin/env python
# -*- coding: utf-8 -*-	

# Filter used by doxygen to any fortran file to do:
# 0 - no error, always return something
# 1 - Add include file (recursive call) 
#     assume that the include file is in the same directory.

import re
import os.path
import sys

#Include command
re_include = re.compile('^[ ]*include.*?\n',re.MULTILINE+re.IGNORECASE)
#Extract file of the include file
re_include_file = re.compile("""^[ ]*include[ ]+['"](?P<file>[\w./_-]+)['"]""",re.IGNORECASE)

def read_file(filename):
    """Read a file and check if an include statement is present"""
    dirname = os.path.dirname(filename)
    if dirname == "": dirname = "."
    for line in open(filename).readlines():
        if re_include.match(line):
            #Name of the include file
            name = re_include_file.match(line).groupdict()['file']
            #Recursive call (very safe)
            try:
                sys.stdout.write("!-- Add include file : %s" % line)
                read_file('%s/%s' % (dirname,name))
            except:
                sys.stdout.write("!-- Failed to include the file %s/%s" % (dirname,name))
                sys.stdout.write(line)
                pass
        else:
            sys.stdout.write(line)


if __name__ == '__main__':
    try:
        read_file(sys.argv[1])
    except:
        pass
    sys.exit(0)
