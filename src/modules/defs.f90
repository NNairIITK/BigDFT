!> @file
!!  File defining parameters for BigDFT package (processed by the build system)
!! @author
!!    Copyright (C) 2008-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!! @warning
!!   THIS FILE IS PROCESSED BY THE BUILD SYSTEM.

!> Modules which contains the low level definitions, as well as some profiling procedures
module module_defs
  use wrapper_MPI
  use wrapper_linalg

  implicit none  

  include 'configure.inc' !< Include variables set from configure.

  integer :: verbose=2    !< Verbosity of the output, control the level of writing (minimal by default)

  ! General precision, density and the wavefunctions types
  integer, parameter :: gp=kind(1.0d0)  !< general-type precision
  integer, parameter :: dp=kind(1.0d0)  !< density-type precision
  integer, parameter :: wp=kind(1.0d0)  !< wavefunction-type precision
  integer, parameter :: tp=kind(1.0d0)  !< DIIS precision (single in this context, if double is only for non-regression)

  !> Define type of data for MPI
  integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypd=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
  !integer, parameter :: mpidtypw=MPI_REAL,mpidtypd=MPI_REAL !in case of single precision

  !> Flag for GPU computing, if CUDA libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical :: GPUconv=.false.,GPUshare=.true.

  !> Flag for GPU computing, if OpenCL libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical :: OCLconv=.false.
  logical :: ASYNCconv=.true.

  !> Logical parameter for the projectors application strategy (true for distributed way)
  !! if the projector allocation passes the memorylimit this is switched to true
  !! inside localize_projectors routines
  logical :: DistProjApply=.true.

  !> experimental variables to test the add of new functionalities
  logical :: experimental_modulebase_var_onlyfion=.false.

  type(mpi_environment) :: bigdft_mpi !< Contains all data needed for MPI processes

  !> Physical constants.
  real(gp), parameter :: Bohr_Ang = 0.5291772108_gp                     !< 1 AU in angstroem
  real(gp), parameter :: Ha_cmm1=219474.6313705_gp                      !< 1 Hartree, in cm^-1 (from abinit 5.7.x)
  real(gp), parameter :: Ha_eV=27.21138386_gp                           !< 1 Hartree, in eV
  real(gp), parameter :: Ha_K=315774.65_gp                              !< 1 Hartree, in Kelvin
  real(gp), parameter :: Ha_THz=6579.683920722_gp                       !< 1 Hartree, in THz
  real(gp), parameter :: Ha_J=4.35974394d-18                            !< 1 Hartree, in J
  real(gp), parameter :: e_Cb=1.602176487d-19                           !< minus the electron charge, in Coulomb
  real(gp), parameter :: kb_HaK=8.617343d-5/Ha_eV                       !< Boltzmann constant in Ha/K
  real(gp), parameter :: amu_emass=1.660538782e-27_gp/9.10938215e-31_gp !< 1 atomic mass unit, in electronic mass
  real(gp), parameter :: AU_GPa=29421.010901602753_gp                   !< 1 Ha/Bohr^3 in GPa

  !> Evergreens
  real(dp), parameter :: pi_param=3.141592653589793238462643383279502884197_dp

  !> Code constants.
  !real(gp), parameter :: UNINITIALISED = -123456789._gp

  !interface for uninitialized variable
  interface UNINITIALIZED
     module procedure uninitialized_dbl,uninitialized_int,uninitialized_real,uninitialized_long
  end interface

  interface
     subroutine bigdft_utils_flush(unit)
       integer, intent(in) :: unit
     end subroutine bigdft_utils_flush
  end interface

  contains

    function uninitialized_int(one) 
      implicit none
      integer(kind = 4), intent(in) :: one
      integer(kind = 4) :: uninitialized_int
      integer :: foo
      foo = kind(one)
      uninitialized_int=-123456789
    end function uninitialized_int

    function uninitialized_long(one) 
      implicit none
      integer(kind = 8), intent(in) :: one
      integer(kind = 8) :: uninitialized_long
      integer :: foo
      foo = kind(one)
      uninitialized_long=-123456789
    end function uninitialized_long

    function uninitialized_real(one) 
      implicit none
      real(kind=4), intent(in) :: one
      real(kind=4) :: uninitialized_real
      integer :: foo
      foo = kind(one)
      uninitialized_real=-123456789.e0
    end function uninitialized_real

    function uninitialized_dbl(one) 
      implicit none
      real(kind=8), intent(in) :: one
      real(kind=8) :: uninitialized_dbl
      integer :: foo
      foo = kind(one)
      uninitialized_dbl=-123456789.d0
    end function uninitialized_dbl

    function fnrm_denpot(x,cplex,nfft,nspden,opt_denpot,user_data)
      use m_ab7_mixing
      implicit none
      !Arguments
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, dimension(*), intent(in) :: x
      integer, dimension(:), intent(in) :: user_data
      !Local variables
      integer :: ierr, ie, iproc, npoints, ishift
      double precision :: fnrm_denpot, ar, nrm_local, dnrm2

      ! In case of density, we use nscatterarr.
      if (opt_denpot == AB7_MIXING_DENSITY) then
         call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
         if (ierr /= 0) then
            call MPI_ABORT(bigdft_mpi%mpi_comm, ierr, ie)
         end if
         npoints = cplex * user_data(2 * iproc + 1)
         ishift  =         user_data(2 * iproc + 2)
      else
         npoints = cplex * nfft
         ishift  = 0
      end if

      ! Norm on spin up and spin down
      nrm_local = dnrm2(npoints * min(nspden,2), x(1 + ishift), 1)
      nrm_local = nrm_local ** 2

      if (nspden==4) then
         ! Add the magnetisation
         ar = dnrm2(npoints * 2, x(1 + cplex * nfft * 2 + ishift), 1)
         ar = ar ** 2
         if (opt_denpot == 0) then
            if (cplex == 1) then
               nrm_local = nrm_local + 2.d0 * ar
            else
               nrm_local = nrm_local + ar
            end if
         else
            nrm_local = 0.5d0 * (nrm_local + ar)
         end if
      end if
      
      ! Summarize on processors
      fnrm_denpot = nrm_local
      if (bigdft_mpi%nproc > 1) then
         call MPI_ALLREDUCE(nrm_local, fnrm_denpot, 1, &
              & MPI_DOUBLE_PRECISION, MPI_SUM, bigdft_mpi%mpi_comm, ierr)
         if (ierr /= 0) call MPI_ABORT(bigdft_mpi%mpi_comm, ierr, ie)
      end if
    end function fnrm_denpot

    function fdot_denpot(x,y,cplex,nfft,nspden,opt_denpot,user_data)
      use m_ab7_mixing
      implicit none
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, intent(in) :: x(*), y(*)
      integer, intent(in) :: user_data(:)

      integer :: ierr, ie, iproc, npoints, ishift
      double precision :: fdot_denpot, ar, dot_local, ddot

      ! In case of density, we use nscatterarr.
      if (opt_denpot == AB7_MIXING_DENSITY) then
         call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
         if (ierr /= 0) then
            call MPI_ABORT(bigdft_mpi%mpi_comm, ierr, ie)
         end if
         npoints = cplex * user_data(2 * iproc + 1)
         ishift  =         user_data(2 * iproc + 2)
      else
         npoints = cplex * nfft
         ishift  = 0
      end if

      if (opt_denpot == 0 .or. opt_denpot == 1) then
         ! Norm on spin up and spin down
         dot_local = ddot(npoints * min(nspden,2), x(1 + ishift), 1, y(1 + ishift), 1)

         if (nspden==4) then
            ! Add the magnetisation
            ar = ddot(npoints * 2, x(1 + cplex * nfft * 2 + ishift), 1, &
                 & y(1 + cplex * nfft * 2 + ishift), 1)
            if (opt_denpot == 0) then
               if (cplex == 1) then
                  dot_local = dot_local + 2.d0 * ar
               else
                  dot_local = dot_local + ar
               end if
            else
               dot_local = 0.5d0 * (dot_local + ar)
            end if
         end if
      else
         if(nspden==1)then
            dot_local = ddot(npoints, x(1 + ishift), 1, y(1 + ishift), 1)
         else if(nspden==2)then
            ! This is the spin up contribution
            dot_local = ddot(npoints, x(1 + ishift + nfft), 1, y(1 + ishift), 1)
            ! This is the spin down contribution
            dot_local = dot_local + ddot(npoints, x(1 + ishift ), 1, y(1 + ishift+ nfft), 1)
         else if(nspden==4)then
            !  \rho{\alpha,\beta} V^{\alpha,\beta} =
            !  rho*(V^{11}+V^{22})/2$
            !  + m_x Re(V^{12})- m_y Im{V^{12}}+ m_z(V^{11}-V^{22})/2
            dot_local = 0.5d0 * (ddot(npoints, x(1 + ishift), 1, y(1 + ishift), 1) + &
                 & dot(npoints, x(1 + ishift), 1, y(1 + ishift + nfft), 1))
            dot_local = dot_local + 0.5d0 * ( &
                 & ddot(npoints, x(1 + ishift + 3 * nfft), 1, y(1 + ishift), 1) - &
                 & ddot(npoints, x(1 + ishift + 3 * nfft), 1, y(1 + ishift + nfft), 1))
            dot_local = dot_local + &
                 & ddot(npoints, x(1 + ishift + nfft), 1, y(1 + ishift + 2 * nfft), 1)
            dot_local = dot_local - &
                 & ddot(npoints, x(1 + ishift + 2 * nfft), 1, y(1 + ishift + 3 * nfft), 1)
         end if
      end if
      
      ! Summarize on processors
      fdot_denpot = dot_local
      if (bigdft_mpi%nproc > 1) then
         call MPI_ALLREDUCE(dot_local, fdot_denpot, 1, &
              & MPI_DOUBLE_PRECISION, MPI_SUM, bigdft_mpi%mpi_comm, ierr)
         if (ierr /= 0) call MPI_ABORT(bigdft_mpi%mpi_comm, ierr, ie)
      end if
    end function fdot_denpot

end module module_defs
