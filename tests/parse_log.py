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
#from yaml_hl import *

#import yaml_hl

class Pouet:
  def __init__(self):
    self.input = "report"
    self.output = None
    self.style = "ascii"
    self.config = os.path.join(path,'yaml_hl.cfg')

# print out a python dictionary in yaml syntax
def dict_dump(dict):
  sys.stdout.write(yaml.dump(dict,default_flow_style=False,explicit_start=True))

#function which removes from a set of lines the yaml_fields contained in the to_remove list
def clean_logfile(logfile_lines,to_remove):
  line_rev=logfile_lines #list of the lines of the logfile
  #loop in the reversed from (such as to parse by blocks)
  extra_lines=20 #internal variable to be customized
  line_rev.reverse()
  #clean the log
  cleaned_logfile=[]
  #for line in line_rev: #line_iter:
  while len(line_rev) >0:
    line=line_rev.pop()
    to_print=line
    #check if the line contains interesting information
    for remove_it in to_remove :
      stream_list=[]
      if remove_it+':' in line.split('#')[0]:
         #creates a new Yaml document starting from the line
         #treat the rest of the line following the key to be removed
         header=''.join(line.split(':')[1:])
         header=header.rstrip()+'\n'
         if len(header) > 0 :
            stream_list.append(header)
         #part to be printed, updated
         to_print = line.split(':')[0] + ": <folded> \n"

         #then check when the mapping will end:
         while True:
            #create a stream with extra_lines block           
            for i in range(0,min(extra_lines,len(line_rev))):
               stream_list.append(line_rev.pop())
            #create a stream to be parsed
            stream=''.join(stream_list)
            #then parse the stream until the last valid position has been found
            try:
              for i in yaml.parse(stream,Loader=yaml.CLoader):
                endpos=i.end_mark.index
            except:
              #convert back the valid stream into a list
              #if needed the stream can be loaded into a document
              item_list=stream[:endpos].split('\n')
              #if lengths are different there is no need to add lines
              if len(item_list) != len(stream_list):
                #last line might be shorter, therefore treat it separately
                last_line=item_list.pop()
                #purge the stream
                for item in item_list:
                  stream_list.remove(item+'\n')
                #extract the remaining line which should be compared with the last one  
                first_line=stream_list.pop(0)[len(last_line.rstrip()):]
                #then put the rest in the line to be treated
                to_print.rstrip('\n')
                to_print += first_line+'\n'
                # the item has been found
                break
         stream_list.reverse()
         #put back the unused part in the document
         line_rev.extend(stream_list)
         # mark that the key has been removed
         print 'removed: ',remove_it
  # then print out the line
    cleaned_logfile.append(to_print)

  return cleaned_logfile


# this is a tentative function written to understand the sanity of the
# BigDFT run
def document_analysis(doc):
  #analyse the energy and the forces
  last=doc["Last Iteration"]
  try:
    last=doc["Last Iteration"]
    dict_dump(last["EKS"])
    stdout.write(yaml.dump(last,default_flow_style=False,explicit_start=True))    
  except:
    last={}
    print 'Last iteration absent'
  try:
    forces=doc["Geometry"]["Forces"]
    dict_dump(forces[0])
    #print forces["maxval"]
  except:
    forces={}
    print 'forces not calculated'
    

def parse_arguments():
  parser = optparse.OptionParser("This script is used to extraact some information from a logfile")
  parser.add_option('-d', '--data', dest='data',default=None, #sys.argv[2],
                    help="BigDFT logfile, yaml compliant (check this if only this option is given)", metavar='FILE')
  parser.add_option('-v', '--invert-match', dest='remove',default=None, #sys.argv[2],
                    help="File containing the keys which have to be excluded from the logfile", metavar='FILE')
  parser.add_option('-o', '--output', dest='output', default="/dev/null", #sys.argv[4],
                    help="set the output file (default: /dev/null)", metavar='FILE')
  #Return the parsing
  return parser


if __name__ == "__main__":
  parser = parse_arguments()
  (args, argtmp) = parser.parse_args()


#args=parse_arguments()
#logfile
with open(args.data, "r") as fp:
  logfile_lines = fp.readlines()
#output file
file_out=open(args.output, "w")
#to_remove list
to_remove = yaml.load(open(args.remove, "r").read(), Loader = yaml.CLoader)
    
#obtain the cleaned document
cleaned_logfile = clean_logfile(logfile_lines,to_remove)

#dump it
for line in cleaned_logfile:
  file_out.write(line) 

sys.exit(0)

#do some tests
document = """---
block sequence:
 - BlockEntryToken
 block mapping:
   ? KeyToken
   : ValueToken
 flow sequence: [FlowEntryToken, FlowEntryToken]
 flow mapping: {KeyToken: ValueToken}
 anchors and tags:
 - &A !!int '5'
 - *A
...
 """

#for token in yaml.scan(document):
#     print token

#ddd
#print args.ref,args.data,args.output

datas    = [a for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader)]
#i=0
#for doc in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader):
#  i+=1
#  print 'doc No.',i#,"Last Iteration" in doc.keys()
  #try:
  #  print yaml.dump(doc["Last Iteration"])
  #except:
  #  print "Last Iteration not found"

#datas    = [a for a in yaml.safe_load_all(open(args.data, "r").read())]
#Profile.run('datas    = [a for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader)]')
#gyi
#runfile = open(args.data, "r").read()
sss
try:
    datas    = [a for a in yaml.load_all(runfile, Loader = yaml.Loader)]
    #Profile.run('datas    = [a for a in yaml.load_all(open(args.data, "r").read())]')
#for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader):
#    print 'New Document'
#    document_analysis(a)
  
except:
  datas = []
  reports = open(args.output, "w")
  fatal_error(args,reports)

ndocs=len(datas)

print 'No. of Documents:',ndocs

for run in datas:
  print 'New Document'
  document_analysis(run)



print 'end'

test = run['Atomic positions within the cell (Atomic and Grid Units)']

#value = test[1]["Ag"]["AU"][1]

#ex={}
#ex["value"]=complex(12,13)
#print 'example',value,3+4j,value,(3+4j),complex(7,8),ex["value"]+complex(value),test[1]["Ag"]["AU"][0]+test[1]["Ag"]["AU"][2]

sys.stdout.write(yaml.dump(run["Last Iteration"],default_flow_style=False,explicit_start=True))
#sys.stdout.write(yaml.dump(run["Direct and transposed data repartition"],default_flow_style=False,explicit_start=True))
#sys.stdout.write(yaml.dump(ex,default_flow_style=False,explicit_start=True))

sys.exit(0)


