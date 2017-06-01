#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import math
import sys
import os,copy,optparse
import cProfile

#path of the script
path=os.path.dirname(sys.argv[0])
#yaml_folder = '/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib'
#if yaml_folder not in sys.path:
#  sys.path.insert(0,'/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib')

import yaml

# print out a python dictionary in yaml syntax
def dict_dump(dict):
  sys.stdout.write(yaml.dump(dict,default_flow_style=False,explicit_start=True))

# this is a tentative function written to extract information from the runs
def document_analysis(doc,to_extract):
  analysis=[]
  for quantity in to_extract:
    #print quantity
    #follow the levels indicated to find the quantity
    value=doc
    for key in quantity:
      #as soon as there is a problem the quantity is null
      try:
        value=value[key]
      except:
        value=None
        break
    if value is not None:
      analysis.append(value)
    else:
      #skip the document is one of the value is None
      return None
  return analysis    
    
def parse_arguments():
  parser = optparse.OptionParser("This script is used to extract some information from a logfile")
  parser.add_option('-d', '--data', dest='data',default=None,action="store_true", #sys.argv[2],
                    help="BigDFT logfile, yaml compliant (check this if only this option is given)", metavar='FILE')
  parser.add_option('-v', '--invert-match', dest='remove',default=None, #sys.argv[2],
                    help="File containing the keys which have to be excluded from the logfile", metavar='FILE')
  parser.add_option('-a', '--analyze', dest='analyze',default=None, #sys.argv[2],
                    help="File containing the keys which have to be analysed and extracted to build the quantities", metavar='FILE')
  parser.add_option('-e', '--extract', dest='extract',default=None, #sys.argv[2],
                    help="File containing the keys which have to be extracted to build the quantities", metavar='FILE')
  parser.add_option('-o', '--output', dest='output', default="/dev/null", #sys.argv[4],
                    help="set the output file (default: /dev/null)", metavar='FILE')
  parser.add_option('-t', '--timedata', dest='timedata',default=False,action="store_true",
                    help="BigDFT time.yaml file, a quick report is dumped on screen if this option is given", metavar='FILE')
  parser.add_option('-n', '--name', dest='name',default=None,
                    help="Give a name to the set of the plot represented", metavar='FILE')
  parser.add_option('-p', '--plot', dest='plottype',default='Seconds',
                    help="Decide the default yscale for the plotting", metavar='FILE')
  parser.add_option('-s', '--static', dest='static',default=False,action="store_true",
                    help="Show the plot statically for screenshot use", metavar='FILE')
  parser.add_option('-f', '--fontsize', dest='fontsize',default=15,
                    help="Determine fontsize of the bar chart plot", metavar='FILE')
  parser.add_option('-k', '--remove-key', dest='nokey',default=False,action="store_true",
                    help="Remove the visualisation of the key from the main plot", metavar='FILE')

  
  #Return the parsing
  return parser

def nopc(string):
  pc=string.rstrip("%")
  try:
    fl=float(pc)
  except:
    fl=-1.0
  return fl

if __name__ == "__main__":
  parser = parse_arguments()
  (args, argcl) = parser.parse_args()


#logfile
#check if timedata is given
if args.timedata:
  import pylab
  from futile import Time
  print 'args of time',args.timedata,argcl
  if args.name is not None:
    title=args.name
  else:
    title='Time bar chart'
  #in the case of more than one file to analyse
  #or in the case of more than one yaml document per file
  #just load the bars data script
  
  #load the first yaml document
  bt=Time.TimeData(*argcl,**args.__dict__)
  bt.show()
  print "hosts",bt.hostnames
  if bt.scf is not None:
    bt.bars_data(title=title) #timing["WFN_OPT"]["Classes"])
    
  if bt.scf[0] is not None and False:
    bt.load_unbalancing(bt.scf[0]["Classes"]) #timing["WFN_OPT"]["Classes"])
      
  if bt.routines[0] is not None and False:
    data=Time.dump_timing_level(bt.routines[0]) #dict_routines)
    plt=Time.polar_axis(data)
    plt.show()
  else:
    pylab.show()

if args.data is None:
  print "No input file given, exiting..."
  exit(0)

if args.analyze is not None and args.data:
  from BigDFT import Logfiles as lf
  instructions=lf.get_log(args.analyze)
  print '#',args.data,argcl
  lf.process_logfiles(argcl,instructions)
  exit(0)

if args.data:
  with open(argcl[0], "r") as fp:
    logfile_lines = fp.readlines()

#output file
file_out=open(args.output, "w")
#to_remove list
if args.remove is not None:
  to_remove = yaml.load(open(args.remove, "r").read(), Loader = yaml.CLoader)
else:
  #standard list which removes long items from the logfile
  to_remove= ["Atomic positions within the cell (Atomic and Grid Units)",
              "Atomic Forces (Ha/Bohr)",
              #"Orbitals",
              #"Energies",
              "Properties of atoms in the system"]
  #to_remove=[]

    
#obtain the cleaned document
from futile import Yaml
cleaned_logfile = Yaml.clean_logfile(logfile_lines,to_remove)

print 'Logfile cleaned, writing output file...'
#dump it
for line in cleaned_logfile:
  file_out.write(line) 

file_out.close()
print 'Output file written'
#extract for any document of the cleaned logfile the final energy and positions
if args.extract is not None:
  to_extract = yaml.load(open(args.extract, "r").read(), Loader = yaml.CLoader)
else:
  to_extract=[ ["Geometry","FORCES norm(Ha/Bohr)","maxval"], ["Last Iteration","EKS"]]

datas=yaml.load_all(''.join(cleaned_logfile), Loader = yaml.CLoader)
extracted_result=[]
try:
  for doc in datas:
    doc_res=document_analysis(doc,to_extract)
    print doc_res,to_extract
    if doc_res is not None:
      extracted_result.append(doc_res)
except:
  pass

print "Number of valid documents:",len(extracted_result)
for it in extracted_result:
  print it
    
iterations = range(len(extracted_result))
energies = [en for [f, en] in extracted_result]
energy_min=min(energies)
energies = [en-energy_min for en in energies]
forces = [f for [f, en] in extracted_result]

import matplotlib.pyplot as plt
#plt.plot(iterations, energies, '.-',label='E - min(E)')
plt.plot(energies, forces, '.-',label='Energy')
#plt.yscale('log')
plt.legend(loc='lower left')
plt.show()
  
sys.exit(0)


