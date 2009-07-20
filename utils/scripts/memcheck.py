#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Copyright (C) 2009 BigDFT group (TD)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#----------------------------------------------------------------------------
# Check malloc.prc (verbose format i.e. memdebug == .true. in memory.f90)
# Date: 11/05/2009
#----------------------------------------------------------------------------

import sys

def array_name(name):
    "Build different structures"
    array = []
    ilast = None
    for i in reversed(name.split("%")[1:]):
        if ilast:
            ilast = i + "%" + ilast
        else:
            ilast = i
        array.append(ilast)
    return array

print "Read the file 'malloc.prc':"
try:
    fd = iter(open("malloc.prc","r").readlines())
except IOError:
    sys.stdout.write("The file 'malloc.prc' does not exist!\n")
    sys.exit(1)

#First line not useful
fd.next()
#Initialized dictionary of variables
variables = dict()
total_size = 0
nalloc = 0
nzero = 0
ndealloc = 0
for line in fd:
    a = line.split()
    #Not used
    routine = a[0]
    #Decomment here to use the fullname and derivatives
    #name = a[1]
    #Here use the last name (after the last "%")
    name = a[1].split("%")[-1]
    if name == "routine":
        #Last line
        continue
    size = int(a[2])
    total_size += size
    if size < 0:
        ndealloc += 1
    elif size > 0:
        nalloc += 1
    else:
        nzero += 1
    if name in variables.keys():
        variables[name][0] += size
        variables[name][1].append((routine,size))
    else:
        variables[name] = [size,[(routine,size)]]

#Group first
keys = variables.keys()
for key in keys:
    for var in array_name(key):
        if var in variables.keys():
            #Group
            variables[var][0] += variables[key][0]
            variables[var][1].append(variables[key][1])
            del variables[key]
            break

print "Remaining memory=%d, allocations=%d, deallocations=%d, zero=%d" % \
    (total_size,nalloc,ndealloc,nzero)

#Check if 0
ok=0
for (key,value) in variables.items():
    if value[0] != 0:
        ok += 1
        print key,value
if ok != 0:
    print "There are %d incoherencies between allocations and deallocations." % ok

