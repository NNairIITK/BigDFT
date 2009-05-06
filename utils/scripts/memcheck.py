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
# Date: 06/05/2009
#----------------------------------------------------------------------------

import sys

def array_name(name):
    'Build different structures'
    ll = name.split('%')
    if len(ll) <= 1:
        return [ name ]
    else:
       var = ll.pop()
       ll.reverse()
       array = [ var ]
       for i in ll:
           var = i + "%" + var
           array.append(var)
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
for line in fd:
    a = line.split()
    routine = a[0]
    name = a[1]
    size = int(a[2])
    if name in variables.keys():
        variables[name][0] += size
        variables[name][1].append((routine,size))
    else:
        variables[name] = [int(a[2]),[(routine,size)]]

#Group first
for key in variables.keys():
    array = array_name(key)
    if len(array) == 1:
        continue
    for var in array:
        if var in variables.keys():
            #Group
            variables[var][0] += variables[key][0]
            variables[var][1].append(variables[key][1])
            del variables[key]

#Check if 0
ok=0
for (key,value) in variables.items():
    if value[0] != 0:
        ok += 1
        print key,value
if ok != 0:
    print "There are %d incoherencies between allocations and deallocations." % ok

