#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

#> @file
## Define a program to run tests to be called by the Makefiles
## @author
##    Copyright (C) 2016-2016 BigDFT group
##    This file is distributed under the terms of the
##    GNU General Public License, see ~/COPYING file
##    or http://www.gnu.org/copyleft/gpl.txt .
##    For the list of contributors, see ~/AUTHORS
#we should think on reducing the compatibility of the argparse to optparse
import UniParse

parser=UniParse.UniParser(description='Regression checker',method='argparse')

parser.option('-f', '--fldiff', dest='fldiff', default="/dev/null",  # sys.argv[4],
              help="set the script file performing fldiff (default: /dev/null)", metavar='FILE')
parser.option('-t', '--tols', dest='tols', default="/dev/null",  # sys.argv[4],
              help="set the yaml file containing the tolerances for each run", metavar='FILE')


#parser.option('-r', '--reference', dest='ref', default=None,  # sys.argv[1],
#              help="reference yaml stream", metavar='REFERENCE')
#parser.option('-d', '--data', dest='data', default=None,  # sys.argv[2],
#              help="yaml stream to be compared with reference", metavar='DATA')
#parser.option('-t', '--tolerances', dest='tols', default=None,  # sys.argv[3],
#              help="File of the tolerances used for comparison", metavar='TOLS')
#parser.option('-l', '--label', dest='label', default=None,
#              help="Define the label to be used in the tolerance file to override the default", metavar='LABEL')
#
args = parser.args()
#print 'Arguments'
print args.__dict__

import os

instr=os.environ.get('F_REGTEST_INSTRUCTIONS')
if instr is None: quit()

import yaml
d_instr=yaml.load(instr)

def run_test(runs):
    for r in runs:
        os.system(r)

def get_time(file):
    if os.path.isfile(file):
        return os.path.getmtime(file)
    else:
        return 0.0


fldiff=args.fldiff
tols=args.tols
base='python '+fldiff+' -t '+tols

for test in d_instr:
    label=test.keys()[0]
    specs=test.values()[0]
    binary=specs.get('binary',label)
    output=specs['output']
    report=label+'.report.yaml'
    ref=specs['reference']
    if get_time(binary) > get_time(output): run_test(specs['runs'])
    os.system(base+' --label '+label+' -r '+ref+' -d '+output+' --output '+report)
        

