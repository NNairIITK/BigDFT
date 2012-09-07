#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import math
import sys
import os,copy,optparse

#path of the script
path=os.path.dirname(sys.argv[0])
#yaml_folder = '/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib'
#if yaml_folder not in sys.path:
#  sys.path.insert(0,'/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib')

import yaml
from yaml_hl import *

start_fail = "<fail>" #"\033[0;31m"
start_fail_esc = "\033[0;31m "
start_success = "\033[0;32m "
start_pass = "<pass>"
end = "</end>" #"\033[m"
end_esc = "\033[m "



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
#  if tols is not None:
#    print 'test',data,ref,tols
  if type(ref) == type({}):
#for a floating point the reference is set for all the lower levels    
    if type(tols) == type(1.0e-1):
      neweps=tols
      tols={}
      for key in ref:
        tols[key]=neweps
    ret = compare_map(data, ref, tols, always_fails)
  elif type(ref) == type([]):
    if type(tols) == type(1.0e-1):
      neweps=tols
      tols=[]
      tols.append(neweps)
    ret = compare_seq(data, ref, tols, always_fails)
  else:
    ret = compare_scl(data, ref, tols, always_fails)
  return ret

#sequence comparison routine
def compare_seq(seq, ref, tols, always_fails = False):
  global failed_checks
  if tols is not None:
    for i in range(len(ref)):
      (failed, newtols) = compare(seq[i], ref[i], tols[0], always_fails)
# Add to the tolerance dictionary a failed result      
      if failed:
        if type(newtols)== type({}):
          tols[0].update(newtols)
        elif type(newtols) == type([]):
          tols[0] = newtols   
        else:
          tols[0] = max(newtols,tols[0])
  else:
    tols = []
    if len(ref) == len(seq):
      for i in range(len(ref)):
        if len(tols) == 0:
          (failed, newtols) = compare(seq[i], ref[i], always_fails = always_fails)
          #  add to the tolerance dictionary a failed result      
          if failed:
            tols.append(newtols)
        else:
          (failed, newtols) = compare(seq[i], ref[i], tols[0], always_fails = always_fails)
          if failed:
            tols[0] = newtols   
    else:
      failed_checks+=1
      if len(tols) == 0:
        tols.append("NOT SAME LENGTH")
      else:
        tols[0] = "NOT SAME LENGTH"
  return (len(tols) > 0, tols)

def compare_map(map, ref, tols, always_fails = False):
  global docmiss,docmiss_it
  if tols is None:
    tols = {}
  for key in ref:
    if not(ignore_key(key)):
      if not(key in map):
        docmiss+=1
        docmiss_it.append(key)
        print "WARNING!!",key,"not found",ref[key]
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
# add to the tolerance dictionary a failed result              
      if failed:
        if key in tols:
          if type(newtols)== type({}):
            tols[key].update(newtols)
          elif type(newtols) == type([]):
            tols[key]=newtols
          else:
            tols[key] = max(newtols,tols[key])
        else:
            tols[key] = newtols
  return (len(tols) > 0, tols)  
  

#compare the scalars.
#return the tolerance if the results are ok
def compare_scl(scl, ref, tols, always_fails = False):
  global failed_checks,discrepancy,biggest_tol
  failed = always_fails
  ret = (failed, None)
#  print scl,ref,tols
#eliminate the character variables
  if type(ref) == type(""):
    if not(scl == ref):
      ret = (True, scl)
  elif not(always_fails):
    if tols is None:
      failed = not(math.fabs(scl - ref) <= epsilon)
    else:
      failed = not(math.fabs(scl - ref) <= tols) 
    discrepancy=max(discrepancy,math.fabs(scl - ref))
    if not(failed):
      if tols is None:
        ret = (always_fails, None)
      else:
        ret = (True, tols)
    else:
      ret = (True, math.fabs(scl - ref))
      if tols is not None:
        biggest_tol=max(biggest_tol,math.fabs(tols))
  if failed:
    failed_checks +=1
  return ret

def document_report(tol,biggest_disc,nchecks,leaks,nmiss,miss_it,timet):

  results={}
  failure_reason = None 

#  disc=biggest_disc
  if nchecks > 0 or leaks > 0 or nmiss > 0:
    if leaks > 0:
      failure_reason="Memory"
    elif nmiss > 0:
      failure_reason="Information"
    else:
      failure_reason="Difference"
  else:
    start = start_success
    message = "succeeded "
    
  results["Test succeeded"]=nchecks == 0  and nmiss==0 and leaks==0
  if failure_reason is not None:
    results["Failure reason"]=failure_reason
  results["Maximum discrepancy"]=biggest_disc
  results["Maximum tolerance applied"]=tol
  results["Seconds needed for the test"]=timet
  if (nmiss > 0):
    results["Missed Reference Items"]=miss_it

  return results


import yaml_hl

class Pouet:
  def __init__(self):
    self.input = "report"
    self.output = None
    self.style = "ascii"
    self.config = os.path.join(path,'yaml_hl.cfg')

def parse_arguments():
  parser = optparse.OptionParser("This script is used to compare two yaml outputs with respect to a tolerance file given. usage: fldiff_yaml.py <options>")
  parser.add_option('-r', '--reference', dest='ref', default=None, #sys.argv[1],
                    help="reference yaml stream", metavar='REFERENCE')
  parser.add_option('-d', '--data', dest='data',default=None, #sys.argv[2],
                    help="yaml stream to be compared with reference", metavar='DATA')
  parser.add_option('-t', '--tolerances', dest='tols', default=None, #sys.argv[3],
                    help="File of the tolerances used for comparison", metavar='TOLS')
  parser.add_option('-o', '--output', dest='output', default="/dev/null", #sys.argv[4],
                    help="set the output file (default: /dev/null)", metavar='FILE')
  parser.add_option('-l', '--label', dest='label', default=None, 
                    help="Define the label to be used in the tolerance file to override the default", metavar='LABEL')

  return parser

if __name__ == "__main__":
  parser = parse_arguments()
  (args, argtmp) = parser.parse_args()


#args=parse_arguments()

#print args.ref,args.data,args.output

references = [a for a in yaml.load_all(open(args.ref, "r"), Loader = yaml.CLoader)]
try:
  datas    = [a for a in yaml.load_all(open(args.data, "r"), Loader = yaml.CLoader)]
except:
  datas = []
  sys.exit(0)
  
orig_tols  = yaml.load(open(args.tols, "r"), Loader = yaml.CLoader)

# take default value for the tolerances
try:
  def_tols = orig_tols["Default tolerances"]
  try:
    keys_to_ignore = orig_tols["Keys to ignore"]
  except:
    keys_to_ignore = []
  try:
    patterns_to_ignore = orig_tols["Patterns to ignore"]
  except:
    patterns_to_ignore = []
except:
  def_tols = {}
  keys_to_ignore = []
  patterns_to_ignore = []

#override the default tolerances with the values which are written in the label
extra_tols={}
if args.label is not None and args.label is not '':
  try:
    extra_tols=orig_tols[args.label]
#adding new keys to ignore    
    try:
      keys_to_ignore += extra_tols["Keys to ignore"]
      del extra_tols["Keys to ignore"]
    except:
      print 'Label',args.label,': No new keys to ignore' 
#adding new patterns to ignore
    try:
      patterns_to_ignore += extra_tols["Patterns to ignore"]
      del extra_tols["Patterns to ignore"]
    except:
      print 'Label',args.label,': No new patterns to ignore'
#adding new tolearnces and override default ones      
    try:
      def_tols.update(extra_tols)
    except:
      print 'Label',args.label,': No new tolerances'
#eliminate particular case  
    del orig_tols[args.label]
  except:
    print 'Label',args.label,' not found in tolerance file'

#determine generic tolerance
try:
  epsilon = orig_tols["Default tolerances"]["Epsilon"]
except:
  epsilon = 0
    
poplist=[]
for k in range(len(keys_to_ignore)):
  if keys_to_ignore[k].find('*') > -1:
    poplist.append(k)
    patterns_to_ignore.append(keys_to_ignore[k].partition('*')[0])
#print poplist
ipopped=0
for k in poplist:
  keys_to_ignore.pop(k-ipopped)
  ipopped+=1

if len(patterns_to_ignore) > 0:
  orig_tols["Patterns to ignore"]=patterns_to_ignore

#print 'Epsilon tolerance',epsilon
#print 'Ignore',keys_to_ignore,'Patterns',patterns_to_ignore

options = Pouet()
failed_documents=0
reports = open(args.output, "w")
max_discrepancy=0.
leak_memory = 0
total_misses=0
total_missed_items=[]
time = 0.
biggest_tol=epsilon
for i in range(len(references)):
  tols={}  #copy.deepcopy(orig_tols)
#  print data
  failed_checks=0
  docmiss=0
  docmiss_it=[]
  discrepancy=0.
  data = datas[i]
  reference = references[i]
#this executes the fldiff procedure
  compare(data, reference, tols)
  try:
    doctime = data["Timings for root process"]["Elapsed time (s)"]
  except:
    doctime = 0
  try:
    docleaks = data["Memory Consumption Report"]["Remaining Memory (B)"]
  except:
    docleaks = 0
  sys.stdout.write("#Document: %2d, failed_checks: %d, Max. Diff. %10.2e, missed_items: %d memory_leaks (B): %d, Elapsed Time (s): %7.2f\n" %\
                  (i, failed_checks,discrepancy,docmiss,docleaks,doctime))
#  print "failed checks",failed_checks,"max diff",discrepancy
  max_discrepancy=max(discrepancy,max_discrepancy)
  #print total time
  time += doctime
  #print remaining memory
  leak_memory += docleaks
  total_misses +=docmiss
  total_missed_items.append(docmiss_it)
  if failed_checks > 0 or docleaks > 0:
    failed_documents+=1
    #optional
    sys.stdout.write(yaml.dump(tols,default_flow_style=False,explicit_start=True))
  newreport = open("report", "w")
  newreport.write(yaml.dump(document_report(biggest_tol,discrepancy,failed_checks,docleaks,docmiss,docmiss_it,doctime),\
                            default_flow_style=False,explicit_start=True))
  newreport.close()
  reports.write(open("report", "rb").read())
  hl = YAMLHighlight(options)
  hl.highlight()
  
#create dictionary for the final report

finres=document_report(biggest_tol,max_discrepancy,failed_documents,leak_memory,total_misses,total_missed_items,time)
if len(references)> 1:
  sys.stdout.write(yaml.dump(finres,default_flow_style=False,explicit_start=True))
  reports.write(yaml.dump(finres,default_flow_style=False,explicit_start=True))

sys.exit(0)


