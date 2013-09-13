from ..importer import modules

BigDFT = modules['BigDFT']

import numpy, ctypes, struct

def __Atoms_str__(at):
  """Defines stringification of an Atoms object."""
  return """<%s.%s at %s (format: %s)>
  - geocode: %c
  - nat    : %d
  - ntypes : %d %s
  - iatype : %s
  - units  : %s
  - alat   : %s
  - coords :\n%s""" % (at.__class__.__module__,
  at.__class__.__name__, hex(id(at)), at.format,
  at.geocode, at.nat, at.ntypes, at.atomnames,
  at.iatype, at.units, at.alat, at.rxyz)

def __Atoms_get_iatype__(at):
  if at.nat == 0:
    return numpy.array([])
  for field_info in at.__class__.__info__.get_fields():
    name = field_info.get_name().replace('-', '_')
    if name == "iatype":
      return numpy.ctypeslib.as_array(ctypes.cast(field_info.get_value(at),
                                                  ctypes.POINTER(ctypes.c_int)),
                                      shape = (at.nat, ))

def __Atoms_get_rxyz__(at):
  if at.nat == 0:
    return numpy.array([])
  for field_info in at.__class__.__info__.get_fields():
    name = field_info.get_name().replace('-', '_')
    if name == "rxyz":
      s = struct.pack('d', field_info.get_value(at))
      return numpy.ctypeslib.as_array(ctypes.cast(struct.unpack('q', s)[0],
                                                  ctypes.POINTER(ctypes.c_double)),
                                      shape = (at.nat, 3))

setattr(BigDFT.Atoms, "__str__", __Atoms_str__)

setattr(BigDFT.Atoms, "iatype", property(fget = __Atoms_get_iatype__))
setattr(BigDFT.Atoms, "rxyz", property(fget = __Atoms_get_rxyz__))

def __Goutput_str__(en):
  """Defines stringification of an Goutput object."""
  return """<%s.%s at %s (etot: %17.17gHa)>
  - Hartree: %14.11gHa
  - XC     : %14.11gHa
  - Vxc    : %14.11gHa
  - ionic  : %14.11gHa
  - disp   : %14.11gHa
  - kinetic: %14.11gHa
  - potent.: %14.11gHa
  - proj.  : %14.11gHa
  - exactxc: %14.11gHa
  - band st: %14.11gHa
  - Kohn-S.: %14.11gHa
  - traceH : %14.11gHa
  - sum(V) : %14.11gHa
  - Vsic   : %14.11gHa
  - pressure: %g Ha/Bohr^3
  - stress  : (Ha/Bohr^3)\n%s
  - forces  : (Ha/Bohr)\n%s""" % (en.__class__.__module__,
  en.__class__.__name__, hex(id(en)), en.etot,
  en.eh, en.exc, en.evxc, en.eion, en.edisp, en.ekin, en.epot, en.eproj, en.eexctX,
  en.ebs, en.eKS, en.trH, en.evsum, en.evsic, en.pressure, en.strten, en.fxyz)

def __Goutput_get_fxyz__(en):
  if en.fdim == 0:
    return numpy.array([])
  for field_info in en.__class__.__info__.get_fields():
    name = field_info.get_name().replace('-', '_')
    if name == "fxyz":
      s = struct.pack('d', field_info.get_value(en))
      return numpy.ctypeslib.as_array(ctypes.cast(struct.unpack('q', s)[0],
                                                  ctypes.POINTER(ctypes.c_double)),
                                      shape = (en.fdim, 3))
    
setattr(BigDFT.Goutput, "__str__", __Goutput_str__)
setattr(BigDFT.Goutput, "fxyz", property(fget = __Goutput_get_fxyz__))
