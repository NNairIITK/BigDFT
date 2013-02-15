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
      s = struct.pack('d', field_info.get_value(at).data)
      return numpy.ctypeslib.as_array(ctypes.cast(struct.unpack('q', s)[0],
                                                  ctypes.POINTER(ctypes.c_double)),
                                      shape = (at.nat, 3))

setattr(BigDFT.Atoms, "__str__", __Atoms_str__)

setattr(BigDFT.Atoms, "iatype", property(fget = __Atoms_get_iatype__))
setattr(BigDFT.Atoms, "rxyz", property(fget = __Atoms_get_rxyz__))


