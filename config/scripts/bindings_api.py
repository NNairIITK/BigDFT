#!/usr/bin/env python

import sys, re

fLine = re.compile(r"^ *(?P<follower>[&]?)(?P<command>[^&!]*) *(?P<continuation>[&]?) *!?.*\n$")

header = None
decl = False
funcLine = int(sys.argv[2].split(":")[1])
iLine = 0
for line in open(sys.argv[2].split(":")[0], "r"):
  iLine += 1
  if iLine == funcLine:
    header = ""
  if header is not None:
    match = fLine.match(line)
    if decl and len(match.group("command")) > 5 and "::" not in match.group("command") and not match.group("command").startswith('!') and len(match.group("follower")) == 0:
      break
    if not match.group("command").startswith("!"):
      header += match.group("command")
      if len(match.group("continuation")) == 0:
        header += '\n'
    if "::" in match.group("command"):
      decl = True

print "/* Fortran header:"
print header.strip()
print "*/"

decl = [x.strip() for x in re.match(r"^[^(]*\((?P<vars>[,a-zA-Z0-9_ ]+)\)", header).group("vars").split(",")]
arg_tp = {}
for line in header.splitlines()[1:]:
  var = line.split("::")
  if not len(var) == 2:
    continue
  tp = [x.strip() for x in var[0].split(",")]
  np = [x.count("(") - x.count(")") for x in tp]
  tp2 = [""]
  acc = 0
  for (ele, n) in zip(tp, np):
    tp2[-1] += ele
    acc += n
    if acc == 0:
      tp2.append("")
    else:
      tp2[-1] += ", "
  tp = tp2
  dim = ""
  for attr in tp:
    if attr.startswith("dimension"):
      ln = len(re.match(r"^[^(]*\((?P<dims>[^)]+)\)", attr).group("dims").split(","))
      if ln > 1:
        dim = "_%dD" % ln
      break
  c_mod = ""
  if tp[0].startswith("type"):
    f_tp = re.match(r"^type\((?P<type>[^)]+)\)", tp[0]).group("type")
    c_tp = "_%s *" % f_tp # "gpointer "
  elif (tp[0].startswith("real") or tp[0].startswith("double")) and "pointer" not in tp:
    c_tp = "double *"
  elif (tp[0].startswith("real") or tp[0].startswith("double")) and "pointer" in tp:
    c_tp = "f90_pointer_double%s " % dim
  elif tp[0].startswith("logical") and "pointer" not in tp:
    c_tp = "int *"
  elif tp[0].startswith("integer") and "pointer" not in tp:
    kind = re.match(r"^[^(]+(\( *[kK][iI][nN][dD] *= *(?P<kind>[48]+) *\))?", tp[0]).group("kind")
    if kind is None or kind == "4":
      c_tp = "int *"
    elif kind == "8":
      c_tp = "long *"
  elif tp[0].startswith("integer") and "pointer" not in tp:
    c_tp = "int *"
  elif tp[0].startswith("integer") and "pointer" in tp:
    c_tp = "f90_pointer_int%s " % dim
  elif tp[0].startswith("character") and "pointer" not in tp:
    c_tp = "char *"
  if "pointer" in tp:
    c_tp += "*"
  if "intent(in)" in tp:
    c_mod = "const "
  args = [re.match(r"^ *(?P<var>[^(]+)", x.strip()).group("var") for x in var[1].split(",")]
  for var in decl:
    if var in args:
      arg_tp[var] = (c_mod, c_tp)

header = ""
if "_" in sys.argv[1]:
  header += "void FC_FUNC_(%s, %s)(" % (sys.argv[1], sys.argv[1].upper())
else:
  header += "void FC_FUNC(%s, %s)(" % (sys.argv[1], sys.argv[1].upper())
ln = len(header)
extra = ""
iExtra = 0
for var in decl:
  (c_mod, c_tp) = arg_tp[var]
  ## if c_tp.startswith("int") and var[0] in "ijklmn":
  ##   header += c_mod + "unsigned " + c_tp
  ## else:
  ##   header += c_mod + c_tp
  header += c_mod + c_tp
  header += var
  header += ", \n"
  header += " " * ln
  if c_tp.startswith("char"):
    iExtra += 1
    extra += ", \n"
    extra += " " * ln
    extra += "int str_ln_%d" % iExtra
header = header[:-3 - ln] + extra + ");"
print header
