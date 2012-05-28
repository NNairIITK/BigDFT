#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import math
import sys, yaml,yaml_hl
import os

start_fail = "<fail>" #"\033[0;31m"
start_success = " " #"<succ>" #"\033[0;32m"
start_pass = "<pass>"
end = "</end>" #"\033[m"


references = [a for a in yaml.load_all(open(sys.argv[1], "r"), Loader = yaml.CLoader)]
datas      = [a for a in yaml.load_all(open(sys.argv[2], "r"), Loader = yaml.CLoader)]
tols       = yaml.load(open(sys.argv[3], "r"), Loader = yaml.CLoader)

#generic document comparison routine
#descend recursively in the dictionary until a scalar is found
#a tolerance value might be passed
def compare(data, ref, tols = None):
  ret = None
  if type(ref) == type({}):
    compare_map(data, ref, tols)
  elif type(ref) == type([]):
    compare_seq(data, ref, tols)
  else:
    ret = compare_scl(data, ref, tols)
  return ret

#sequence comparison routine
def compare_seq(seq, ref, tols):
  if tols is not None and len(tols) == len(ref):
    for i in range(len(ref)):
      ret = compare(seq[i], ref[i], tols[i])
      if ret is not None:
        ref[i] = ret
  else:
    for i in range(len(ref)):
      ret = compare(seq[i], ref[i])
      if ret is not None:
        ref[i] = ret

def compare_map(map, ref, tols):
  for (key, value) in ref.items():
#    print "#####", key
    if type(tols) == type({}) and key in tols:
      ret = compare(map[key], value, tols[key])
    else:
      ret = compare(map[key], value)
    if ret is not None:
      ref[key] = ret

#compare the scalars.
#return the tolerance if the results are ok
def compare_scl(scl, ref, tols):
  global failed_checks
  if type(ref) == type(""):
    if scl == ref:
      ret = ref
    else:
      failed_checks +=1
      ret = start_fail + scl + end
  else:
    if tols is None:
      ret = math.fabs(scl - ref) == 0
    else:
      ret = math.fabs(scl - ref) < tols
    if ret:
      if tols is None:
        ret = ".inf" #"%s.inf%s" % (start_success, end) #float("inf") #
      else:
        ret = tols
    else:
      failed_checks +=1
      ret = start_fail + str(scl - ref)  + " (" + str(tols) + ")" + end
#  print ret
  return ret


failed_documents=0
for i in range(len(references)):
#  print data
  failed_checks=0
  data = datas[i]
  reference = references[i]
  print "Document ",i,":"
  #sys.stdout.write(yaml.dump(data).encode('utf8'))
  compare(data, reference, tols)
  print "failed checks",failed_checks
  if failed_checks > 0:
    failed_documents+=1
    sys.stdout.write(yaml.dump(reference))
#  print yaml.scan(reference)
#  print substitute(reference,tols)
#  hl.input=reference
#  hl.output=stdout
if failed_documents > 0:
  sys.stderr.write("Test failed for %s Documents\n" % failed_documents)
  sys.exit(1)
else:
  print 'Test passeed'
  sys.exit(0)


  #print reference
#  sys.stdout.write(yaml.dump(reference).encode('utf8'))

