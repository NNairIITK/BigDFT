#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import sys, yaml
import os

fin = open(os.path.dirname(sys.argv[1]) + "/data/time.yaml", "r")

for data in yaml.load_all(fin):
  print data

fin.close()
