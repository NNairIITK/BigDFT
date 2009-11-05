#!/usr/bin/env python
# -*- coding: us-ascii -*-

import re
import os
import sys

#Values to write in green
start = "\033[0;32m"
end = "\033[m"

def do_cache(ncache,args):
    "Do a calculation to test a value of ncache"
    text = "integer, parameter :: ncache"
    re_ncache = re.compile(text+".*",re.MULTILINE)
    file_module = open("fft3d-temp.f90").read()
    file_module = re_ncache.sub(text+"=%d" % ncache,file_module)
    open("fft3d-temp.f90","w").write(file_module)
    sys.stdout.write("compile.")
    sys.stdout.flush()
    os.system("%s -c fft3d-temp.f90" % args)
    os.system("%s -o fft_cache fft_cache-temp.f90 fft3d-temp.f90" % args)
    sys.stdout.write("test")
    sys.stdout.flush()
    #os.system("./fft_cache >> fft_cache.out 2>&1")
    os.system("./fft_cache")

#List of tested values of ncache
list_cache = [ 0, 6, 12, 24, 50, 75, 100]
#list_cache = range(0,2000,10)

if len(sys.argv) == 1:
    sys.stderr.write("Usage: fft_cache name_of_compiler_and_options\n")
    sys.exit(1)
args=""
for i in sys.argv[1:]:
    args += "%s " % i

#Remove the file "fft_cache.dat"
if os.path.exists("fft_cache.dat"):
    os.remove("fft_cache.dat")
if os.path.exists("fft_cache.out"):
    os.remove("fft_cache.out")

for ncache in list_cache:
    sys.stdout.write("[%d:" % ncache)
    do_cache(ncache*1024,args)
    sys.stdout.write("]")
sys.stdout.write("done.\n")

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

