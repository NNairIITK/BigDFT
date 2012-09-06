#!/usr/bin/env python

import os
import re

codes = {1:"1", 11:"2"}

def export_psp(dir, symbol, elements):
  f = open(os.path.join(dir, symbol), "r")
  f.readline()
  psppar = []
  (nzatom, nelpsp) = map(int, f.readline().split()[0:2])
  (pspcod, ixcpsp) = map(int, f.readline().split()[0:2])
  # read(11,*) psppar(0,0),nn,(psppar(0,j),j=1,nn) !local PSP parameters
  vals = f.readline().split()
  psppar.append([float(vals[0])] + map(float, vals[2:2 + int(vals[1])]))
  # read(11,*) nlterms !number of channels of the pseudo
  for l in range(int(f.readline().split()[0])):
    # h_ij terms
    vals = f.readline().split()
    size = int(vals[1])
    coeffs = map(float, vals[2:2 + size]) + [0.] * (3 - size)
    for j in range(size - 1):
      coeffs += map(float, f.readline().split()[:size - j - 1]) + [0.] * (3 - size)
    for j in range(size, 3):
      coeffs += [0.] * (3 - j)
    hij = [coeffs[0], coeffs[3], coeffs[5], coeffs[1], coeffs[2], coeffs[4]]
    # k_ij terms
    if l > 0:
      coeffs = map(float, f.readline().split()[:size]) + [0.] * (3 - size)
      for j in range(size - 1):
        coeffs += map(float, f.readline().split()[:size - j - 1]) + [0.] * (3 - size)
      for j in range(size, 3):
        coeffs += [0.] * (3 - j)
      kij = [coeffs[0], coeffs[3], coeffs[5], coeffs[1], coeffs[2], coeffs[4]]
    psppar.append([float(vals[0])] + hij)

  f.close()

  ele = symbol.split("-")[0]
  if ele in elements:
    ele += "_sc"    
  if ele in elements:
    ele += "+"
  elements.add(ele)
  print '  else if (trim(symbol) == "%s" .and. trim(name_ixc) == trim(name_xcpsp(%s))) then' % (ele, codes[ixcpsp])
  print '     nzatom   = %d' % nzatom
  print '     nelpsp   = %d' % nelpsp
  print '     npspcode = %d' % pspcod
  print '     ixc      = %d' % ixcpsp
  l = 0
  for vals in psppar:
    print '     psppar(%d,0:%d) = (/ %12.12f_gp' % (l, len(vals) - 1, vals[0]),
    k = 0
    for val in vals[1:]:
      print ', %12.12f_gp' % val,
      if (k == 3):
        print ' &\n        & ',
      k += 1
    print ' /)'
    l += 1
  print '     exists = .true.'

def natural_sortkey(string):
  ret = []
  for num, alpha in tokenize(string):
    if num:
      ret.append(int(num))
    else:
      ret.append(alpha)
  return ret
## Disabled since not working for Python 2.4 ...
##  return tuple(int(num) if num else alpha for num, alpha in tokenize(string))

re_psp = re.compile("[A-Z][a-z]?-q[0-9]+")
tokenize = re.compile(r'(\d+)|(\D+)').findall

elements = set()
pspdir = "utils/PSPfiles/Krach-LDA"
for f in sorted(os.listdir(pspdir), key = natural_sortkey):
  if re_psp.match(f) is not None:
    export_psp(pspdir, f, elements)

elements = set()
pspdir = "utils/PSPfiles/Krach-PBE"
for f in sorted(os.listdir(pspdir), key = natural_sortkey):
  if re_psp.match(f) is not None:
    export_psp(pspdir, f, elements)
