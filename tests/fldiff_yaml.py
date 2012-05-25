#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import sys, yaml

fin = open(sys.argv[1], "r")

pin = yaml.parse(fin)
for e in pin:
  print e

fin.close()
