#!/usr/bin/env python
# -*- coding: us-ascii -*-

import re
import os
import sys

#Values to write in green
start = "\033[0;32m"
end = "\033[m"

#List of tested values of ncache
list_cache = [ 0, 6, 12, 24, 50, 75, 100]
#list_cache = range(0,2000,10)

#Build new file
dico = dict()
caches = set()
lengthes = set()
for line in open("fft_cache.dat").readlines():
    line = line.split()
    if len(line) > 2:
        cache = int(line[0])
        caches.add(cache)
        length = int(line[1])
        lengthes.add(length)
        if not dico.has_key(length):
            dico[length] = dict()
        dico[length][cache] = float(line[-1])

caches = list(caches)
caches.sort()
lengthes = list(lengthes)
lengthes.sort()

fd = open("fft_columns.dat","w")
fd.write("#n1")
for c in caches:
    fd.write(" %d" % c)
fd.write("\n")

performances = dict()
for c in caches:
    performances[c] = 0.

for l in lengthes:
    fd.write("%d " % l)
    best = 1.e10
    for c in caches:
        if dico[l].has_key(c):
            perf = dico[l][c]
            fd.write("%f20.15 " % perf)
            if perf < best:
                best = perf
            if l >= 50 and l <= 500:
                performances[c] += perf
        else:
            fd.write(" * ")
    fd.write("%f20.15\n" % best)
fd.close()

best = caches[0]
fd = open("fft_perf.dat","w")
fd.write("#ncache performances (50 <= x <+ 500) time\n")
for c in caches:
    perf = performances[c] 
    if performances[c] < performances[best]:
        best = c
    fd.write("%s %s\n" % (c,perf))
fd.close()

print "Use fft_cache.gnuplot to display the results"
print "The best ncache between tested values for 50 <= n1 <= 500 is ncache=",best
print start+"Test succeeded"+end

