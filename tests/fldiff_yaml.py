#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import math
import sys, yaml
import os
from yaml_hl import *

start_fail = "<fail>" #"\033[0;31m"
start_success = " " #"<succ>" #"\033[0;32m"
start_pass = "<pass>"
end = "</end>" #"\033[m"


references = [a for a in yaml.load_all(open(sys.argv[1], "r"), Loader = yaml.CLoader)]
datas      = [a for a in yaml.load_all(open(sys.argv[2], "r"), Loader = yaml.CLoader)]
tols       = yaml.load(open(sys.argv[3], "r"), Loader = yaml.CLoader)

# take default value for the tolerances
try:
  def_tols = tols["Default tolerances"]
  try:
    epsilon = tols["Default tolerances"]["Epsilon"]
  except:
    epsilon = 0
  try:
    keys_to_ignore = tols["Keys to be ignored"]
  except:
    keys_to_ignore = []
except:
  def_tols = {}
  epsilon = 0
  keys_to_ignore = []

patterns_to_ignore = []

poplist=[]
for k in range(len(keys_to_ignore)):
  if keys_to_ignore[k].find('*') > -1:
    poplist.append(k)
    patterns_to_ignore.append(keys_to_ignore[k].partition('*')[0])
print poplist
for k in poplist:
  keys_to_ignore.pop(k)

print 'Ignore',keys_to_ignore,'Patterns',patterns_to_ignore

def ignore_key(key):
  ret = key in keys_to_ignore
  if (not(ret)):
    for p in patterns_to_ignore:
      if key.find(p) > -1:
        ret=True
        exit
  return ret
    

#generic document comparison routine
#descend recursively in the dictionary until a scalar is found
#a tolerance value might be passed
def compare(data, ref, tols = None, always_fails = False):
  if type(ref) == type({}):
    ret = compare_map(data, ref, tols, always_fails)
  elif type(ref) == type([]):
    ret = compare_seq(data, ref, tols, always_fails)
  else:
    ret = compare_scl(data, ref, tols, always_fails)
  return ret

#sequence comparison routine
def compare_seq(seq, ref, tols, always_fails = False):
  if tols is not None:
    for i in range(len(ref)):
      (failed, newtols) = compare(seq[i], ref[i], tols[0], always_fails)
      if failed:
        tols[0] = newtols
  else:
    tols = []
    for i in range(len(ref)):
      (failed, newtols) = compare(seq[i], ref[i], always_fails = always_fails)
      if failed:
        if len(tols) == 0:
          tols.append(newtols)
        else:
          tols[0] = newtols   
  return (len(tols) > 0, tols)

def compare_map(map, ref, tols, always_fails = False):
  if tols is None:
    tols = {}
  for key in ref:
#    print 'key',key,ignore_key(key)
    if not(ignore_key(key)):
      if not(key in map):
        print key,"not found",ref[key]
        always_fails = True
        value = ref[key]
      else:
        value = map[key]
      if type(tols) == type({}) and key in tols:
        (failed, newtols) = compare(value, ref[key], tols[key], always_fails)
      elif key in def_tols: #def tols is rigorously a one-level only dictionary
        (failed, newtols) = compare(value, ref[key], def_tols[key], always_fails)
      else:
        (failed, newtols) = compare(value, ref[key], always_fails = always_fails)
      if failed:
#        print 'here,newtols',newtols
        tols[key] = newtols
  return (len(tols) > 0, tols)  
  

#compare the scalars.
#return the tolerance if the results are ok
def compare_scl(scl, ref, tols, always_fails = False):
  global failed_checks,max_discrepancy
  failed = always_fails
  ret = (failed, None)
  if type(ref) == type(""):
    if not(scl == ref):
      ret = (True, scl)
  elif not(always_fails):
    max_discrepancy=max(max_discrepancy,math.fabs(scl - ref))
    if tols is None:
      failed = not(math.fabs(scl - ref) < epsilon)
    else:
#      print math.fabs(scl - ref), float(str(tols).split()[0])
#      print math.fabs(scl - ref) - float(str(tols).split()[0])
      failed = not(math.fabs(scl - ref) <= float(str(tols).split()[0]))
    if not(failed):
      if tols is None:
        ret = (always_fails, None)
      else:
        ret = (True, float(str(tols).split()[0]))
    else:
      ret = (True, str(math.fabs(scl - ref)) + " FAILED")
  if failed:
    failed_checks +=1
#  print ret
  return ret

import yaml_hl

class Pouet:
  def __init__(self):
    self.input = "report"
    self.output = None
    self.style = "ascii"
    self.config = 'yaml_hl.cfg'

options = Pouet()
failed_documents=0
reports = open("reports", "w")
for i in range(len(references)):
#  print data
  failed_checks=0
  max_discrepancy=0.
  data = datas[i]
  reference = references[i]
  #sys.stdout.write(yaml.dump(data).encode('utf8'))
  compare(data, reference, tols)
  sys.stdout.write("#Document: %d, failed_checks: %d\n" % (i, failed_checks))
  print "failed checks",failed_checks,"max diff",max_discrepancy
  if failed_checks > 0:
    newtols = open("report", "w")
    newtols.write("#Document: %d, failed_checks: %d\n" % (i, failed_checks))
    failed_documents+=1
    newtols.write(yaml.dump(tols))
    newtols.close()
    reports.write(open("report", "rb").read())
    hl = YAMLHighlight(options)
    hl.highlight()
  
#  print yaml.scan(reference)
#  print substitute(reference,tols)
#  hl.input=reference
#  hl.output=stdout


if failed_documents > 0:
  sys.stderr.write("Test failed for %s Documents\n" % failed_documents)
  sys.exit(1)
else:
  print 'Test passed.'
  sys.exit(0)


  #print reference
#  sys.stdout.write(yaml.dump(reference).encode('utf8'))

