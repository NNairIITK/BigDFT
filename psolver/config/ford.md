---
Title: PSolver documentation
Summary: How to use the Poisson Solver library.
docmark_alt: '<'
source: false
graph: false:
search: false
exclude: Build_Kernel.f90
exclude: conv_check_fft.f90
exclude: conv_check_fftw.f90
exclude: createKernel.f90
exclude: cufft_fake.f90
exclude: environment.f90
exclude: exctx_calculation.f90
exclude: FDder.f90
exclude: FDder_nehalem.f90
exclude: FDder_nehalem_kernels.f90
exclude: IObox.f90
exclude: IOetsf.f90
exclude: IOetsf_fake.f90
exclude: mp_quadrature.f90
exclude: PSbase.f90
exclude: PSbox.f90
exclude: PSolver_Base_new.f90
exclude: PSolver_Main.f90
exclude: PStypes.f90
exclude: scaling_function.f90
exclude: wofz.f90
---

@note
This Poisson Solver is part of the BigDFT project but can also be used separately.
All the source files of this directory are delivered under the GNU General Public License; 
please see the file COPYING for copying conditions.

See the file INSTALL for generic compilation and installation instructions.

##References
This Poisson Solver code is built following the method described in the papers:

- G. Fisicaro, L. Genovese, O. Aindreussi, N. Marzari, S. Goedecker
    + Journal of Chemical Physics, 144 (1), 014103 (2016)
    + A generalized Poisson and Poisson-Boltzmann solver for electrostatic environments
    + http://dx.doi.org/10.1063/1.4939125
- A. Cerioni, L. Genovese, A. Mirone, V. A. Sole
    + Journal of Chemical Physics, 137, 134108 (2012)
    + Efficient and accurate solver of the three-dimensional screened and unscreened Poisson's equation with generic boundary conditions
    + http://dx.doi.org/10.1063/1.4755349
- L.Genovese, T. Deutsch, S. Goedecker,
    + Journal of Chemical Physics, 127 (5), 054704 (2007).
    + Efficient and accurate three-dimensional Poisson solver for surface problems.
    + http://dx.doi.org/10.1063/1.2754685
- L. Genovese, T. Deutsch, A. Neelov, S. Goedecker, G. Beylkin,
    + Journal of Chemical Physics, 125 (7), 074105 (2006).
    + Efficient solution of Poisson's equation with free boundary conditions
    + http://dx.doi.org/10.1063/1.2335442

Citing of this references is greatly appreciated if the routines are used 
for scientific work.
