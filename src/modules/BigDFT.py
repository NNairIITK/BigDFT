from ..importer import modules
from ..overrides import override

BigDFT = modules['BigDFT']

import numpy, ctypes, struct

__all__ = []

def __is_dict_scalar__(v):
  return isinstance(v, int) or isinstance(v, float) or isinstance(v, str) or isinstance(v, bool)

class DictAccessor:
  def __init__(self, dict, position, elems = slice(0, -1, 1)):
    self.dict = dict
    self.position = position
    self.elems = elems
    
  def __getitem__(self, k):
    self.dict.move_to(self.position)
    if isinstance(k, str):
      (valid, it) = self.dict.move_to_key(k)
    elif isinstance(k, int):
      if k < 0:
        k += self.dict.len()
      (valid, it) = self.dict.move_to_item(k)
    elif isinstance(k, slice):
      return DictAccessor(self.dict, self.position, k)
    else:
      raise TypeError(type(k))
    if not(valid):
      raise IndexError(k)
    return DictAccessor(self.dict, it)

  def __setitem__(self, k, v):
    self.dict.move_to(self.position)
    if isinstance(k, str):
      self.dict.pop(k)
      self.dict.update(v, self.dict.insert(k))
    elif isinstance(k, int):
      (valid, it) = self.dict.move_to_item(k)
      if not(valid):
        raise IndexError(k)
      self.dict.update(v, it)
    else:
      raise AttributeError

  def __iter__(self):
    if not(self.dict.value() == "__list__"):
      raise TypeError("Not iterable")
    self.dict.move_to(self.position)
    self.indices = self.elems.indices(self.dict.len())
    self.id = 0
    return self

  def next(self):
    i = self.indices[0] + self.indices[2] * self.id
    if i > self.indices[1]:
      self.dict.move_to(self.position)
      raise StopIteration
    # Move to root of the list.
    self.dict.move_to(self.position)
    (valid, it) = self.dict.move_to_item(i)
    if not(valid):
      self.dict.move_to(self.position)
      raise StopIteration
    self.id += 1
    return DictAccessor(self.dict, it)

  def map(self, func):
    self.dict.move_to(self.position)
    val = self.dict.value()
    if val == "__list__":
      val = []
      for ele in self:
        val.append(ele.map(func))
    elif val == "__dict__":
      val = {}
      (valid, current) = self.dict.iter()
      while (valid):
        # the func() call may move the current pointer of dict.
        k = self.dict.key()
        val[k] = DictAccessor(self.dict, current).map(func)
        self.dict.move_to(current)
        (valid, current) = self.dict.next()
    else:
      val = func(val)
    return val

  def __str__(self):
    return str(self.map(str))

  def __mul__(self, val):
    if isinstance(val, int):
      return numpy.array(self.map(int)) * val
    elif isinstance(val, float):
      return numpy.array(self.map(float)) * val
    else:
      raise TypeError

  def __len__(self):
    self.dict.move_to(self.position)
    return int(self.dict.len())

class Dict(BigDFT.Dict):
  def __new__(cls, args = (), kwargs = {}):
    (obj, root) = BigDFT.Dict.new()
    obj.rootIter = root
    return obj

  def __init__(self, source = None):
    if source is not None:
      self.update(source, self.rootIter)

  def __dict_add__(self, var, it):
    for (k, v) in var.items():
      self.move_to(it)
      self.update(v, self.insert(k))

  def __list_add__(self, var, it):
    for v in var:
      self.move_to(it)
      self.update(v, self.append())

  def update(self, add, it):
    if __is_dict_scalar__(add):
      # scalar case
      self.set(None, str(add))
    elif isinstance(add, dict):
      # dictionary case
      self.__dict_add__(add, it)
    elif hasattr(add, "__iter__"):
      # List case
      self.__list_add__(add, it)
    else:
      raise TypeError

  def __getitem__(self, k):
    return DictAccessor(self, None)[k]

  def __setitem__(self, k, v):
    DictAccessor(self, None)[k] = v
    
Dict = override(Dict)
__all__.append('Dict')

class Atoms(BigDFT.Atoms):
  def __str__(self):
    """Defines stringification of an Atoms object."""
    return """<%s.%s at %s (format: %s)>:
    - geocode: %c
    - nat    : %d
    - ntypes : %d %s
    - iatype : %s
    - units  : %s
    - alat   : %s
    - coords :\n%s""" % (self.__class__.__module__,
    self.__class__.__name__, hex(id(self)), self.format,
    self.geocode, self.nat, self.ntypes, self.atomnames,
    self.iatype, self.units, self.alat, self.rxyz)

  @property
  def iatype(self):
    if self.nat == 0:
      return numpy.array([])
    return numpy.ctypeslib.as_array(ctypes.cast(super(BigDFT.Atoms, self).iatype,
                                                ctypes.POINTER(ctypes.c_int)),
                                    shape = (self.nat, ))

  @property
  def ifrztyp(self):
    if self.nat == 0:
      return numpy.array([])
    return numpy.ctypeslib.as_array(ctypes.cast(super(BigDFT.Atoms, self).ifrztyp,
                                                ctypes.POINTER(ctypes.c_int)),
                                    shape = (self.nat, ))

  @property
  def rxyz(self):
    if self.nat == 0:
      return numpy.array([])
    s = struct.pack('d', super(BigDFT.Atoms, self).rxyz)
    return numpy.ctypeslib.as_array(ctypes.cast(struct.unpack('q', s)[0],
                                                ctypes.POINTER(ctypes.c_double)),
                                    shape = (self.nat, 3))

  @property
  def format(self):
    return ''.join(map(chr, super(BigDFT.Atoms, self).format))

  @property
  def units(self):
    return ''.join(map(chr, super(BigDFT.Atoms, self).units))

Atoms = override(Atoms)
__all__.append('Atoms')

def __Goutput_str__(en):
  """Defines stringification of an Goutput object."""
  return """<%s.%s at %s (etot: %17.17gHa)>:
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

class Inputs(BigDFT.Inputs):
  def __new__(cls, naming = None, dictionary = {}):
    obj = BigDFT.Inputs.new(naming)
    if len(dictionary) > 0:
      for pair in dictionary.items():
        obj.set(*pair)
    return obj

  def set(self, *args):
    if len(args) == 1 and isinstance(args[0], dict):
      for pair in args[0].items():
        self.set(*pair)
      return
    (key, val) = args
    if isinstance(val, tuple) or isinstance(val, list):
      if all([isinstance(elem, tuple) or isinstance(elem, list) for elem in val]):
        for (i, elem) in enumerate(val):
          self.set_array_at(key, i, map(str, elem))
      else:
        self.set_array(key, map(str, val))
    else:
      super(BigDFT.Inputs, self).set(key, str(val))

Inputs = override(Inputs)
__all__.append('Inputs')

class Run(BigDFT.Run):
  def __new__(cls, source, dump = True):
    if isinstance(source, str):
      var = BigDFT.Dict.new_from_yaml(source)
    elif isinstance(source, dict):
      var = BigDFT.Dict(source)
    else:
      raise TypeError
    return BigDFT.Run.new_from_dict(var, dump)

  def __init__(self, args = (), kwargs = {}):
    return
      
Run = override(Run)
__all__.append('Run')
