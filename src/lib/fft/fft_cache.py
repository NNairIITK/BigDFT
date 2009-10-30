#!/usr/bin/env python
# -*- coding: us-ascii -*-

import re
import os

text = "integer, parameter :: ncache"

def do_cache(ncache):
    re_ncache = re.compile(text+".*",re.MULTILINE)
    file_module = open("fft3d.f90").read()
    file_module = re_ncache.sub(text+"=%d\n" % ncache,file_module)
    open("toto.f90","w").write(file_module)
    os.system("gfortran -c toto.f90")
    os.system("gfortran -o fft_cache fft_cache.f90 toto.f90")
    os.system("./fft_cache")

list = [ 1, 6, 10, 20, 50, 100, 500, 1000]

#Remove the file "fft_cache.dat"
os.remove("fft_cache.dat")
for ncache in list:
    do_cache(ncache*1024)
    print ncache

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
        dico[length][cache] = line[-1]

caches = list(caches)
caches.sort()
lengthes = list(lengthes)
lengthes.sort()

fd = open("fft_columns.dat","w")
fd.write("#n1")
for c in caches:
    fd.write(" %d" % c)
fd.write("\n")

for l in lengthes:
    fd.write("%d" % l)
    for c in caches:
        if dico[l].has_key(c):
            fd.write("%s" % dico[l][c])
        else:
            fd.write(" * ")
    fd.write("\n")
fd.close()
