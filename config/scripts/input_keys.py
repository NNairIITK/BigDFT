#!/usr/bin/env python

import sys

keys = []
for line in open("src/modules/input_keys.f90", "r").xreadlines():
  if line.strip() == "contains":
    break
  if "::" in line:
    decls = line.split("::")[1]
    decls = decls.split(',')
    np = [x.count('"') + x.count("'") for x in decls]
    decls2 = [""]
    acc = 0
    for (decl, n) in zip(decls, np):
      decls2[-1] += decl
      acc += n
      if acc % 2 == 0:
        decls2.append("")
      else:
        decls2[-1] += ","
    decls = decls2[:-1]
    for decl in decls:
      try:
        (name, value) = [x.strip() for x in decl.split("=", 1)]
        keys.append((name, value[1:-1]))
      except:
        continue

strs = {}

strs["static"] = """#ifndef INPUT_KEYS_H
#define INPUT_KEYS_H

static const gchar* _input_keys[] = {
"""
strs["static_add"] = ""
strs["enum"] = """typedef enum {\n"""

for (name, value) in keys[:-1]:
  if name.endswith("_VARIABLES"):
    tp = name
  strs["static"]     += "  \"%s\",\n" % value
  strs["enum"]       += "  INPUTS_%s,\n" % name
  strs["static_add"] += "  INPUTS_%s,\n" % tp
strs["static"]     += "  \"%s\"" % keys[-1][1]
strs["enum"]       += "  INPUTS_%s" % keys[-1][0]
strs["static_add"] += "  INPUTS_%s" % tp

strs["static"] += """
};

static const BigDFT_InputsKeyIds _input_files[] = {
%s
};

#endif""" % strs["static_add"]
strs["enum"] += """
} BigDFT_InputsKeyIds;"""

if len(sys.argv) == 1 or sys.argv[1] == "static":
  print strs["static"]
if len(sys.argv) == 1 or sys.argv[1] == "enum":
  print strs["enum"]
