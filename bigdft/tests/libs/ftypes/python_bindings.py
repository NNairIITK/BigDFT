#!/usr/bin/env python

from gi.repository import BigDFT

import numpy, ctypes

at = BigDFT.Atoms.new_from_file("posinp.ascii")
at.rxyz[1, 2] = 10.
print at

at2 = BigDFT.Atoms.new()
print at2
at2.set_n_types(at.ntypes)
at2.set_n_atoms(at.nat)
at2.iatype[:] = at.iatype[:]
at2.alat[:] = at.alat[:]
at2.rxyz[:,:] = at.rxyz[:,:]
print at2
