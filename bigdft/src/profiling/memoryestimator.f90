!> @file
!!  Routines to estimate the use of memory
!! @author
!!    Copyright (C) Luigi Genovese, CEA Grenoble, France, 2007-2013
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Estimation of the used memory
subroutine MemoryEstimator(nproc,idsx,lr,norb,nspinor,nkpt,nprojel,nspin,itrpmax,iscf,mem)

  use module_base
  use module_types
  use Poisson_Solver
  use locreg_operations, only: memspace_work_arrays_sumrho,memspace_work_arrays_locham
  implicit none

  !Arguments
  integer, intent(in) :: nproc,idsx,norb,nspin,nprojel
  integer, intent(in) :: nkpt,nspinor,itrpmax,iscf
  type(locreg_descriptors), intent(in) :: lr
  type(memory_estimation), intent(out) :: mem
  !Local variables
  !character(len=*), parameter :: subname='MemoryEstimator'
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: norbp,nvctrp,n1,n2,n3
  integer :: n01,n02,n03,m1,m2,m3,md1,md2,md3,nd1,nd2,nd3
  integer(kind=8) :: mworkham, mworkrho
  real(kind=8) :: omemwf,omemker,omemden,omempot,omemproj,nden,npotden,npotham,narr
  real(kind=8) :: tt
!!$ real(kind=8) :: timinamount

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3

  !here we must add the estimation for the projectors

  tt=dble(norb * nkpt)/dble(nproc)
  norbp=int((1.d0-eps_mach*tt) + tt)
  tt=dble(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)/dble(nproc)
  nvctrp=int((1.d0-eps_mach*tt) + tt)

  ! Multiply the size of one wavefunction per its spinor value.
!!$  if(nspin==4) then !quadruple size for complex spinors
  norbp=norbp ! do not multiply also the number of bands.
  nvctrp=nvctrp*nspinor
!!$  end if

  !wavefunction memory per orbitals
  omemwf=real(nvctrp*nproc*8,kind=8)
  
  if (lr%geocode == 'P') then
     call P_FFT_dimensions(2*n1+2,2*n2+2,2*n3+2,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
     n01=2*n1+2
     n02=2*n2+2
     n03=2*n3+2
  else if (lr%geocode == 'S') then
     call S_FFT_dimensions(2*n1+2,2*n2+31,2*n3+2,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
     n01=2*n1+2
     n02=2*n2+31
     n03=2*n3+2
  else if (lr%geocode == 'F') then
     call F_FFT_dimensions(2*n1+31,2*n2+31,2*n3+31,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
     n01=2*n1+31
     n02=2*n2+31
     n03=2*n3+31
  end if
  tt = 8.d0*real(n1*n2*n3,kind=8)/real(n01*n02*n03,kind=8)

  !density memory
  omemden=real(md3*md2/nproc,kind=8)*8.d0*real(md1*nspin,kind=8)
  !kernel memory
  omemker=real(nd2*nd3/nproc,kind=8)*8.d0*real(nd1,kind=8)
  !memory of full grid arrays
  omempot=real(n02*n03,kind=8)*8.d0*real(n01*nspin,kind=8)
  !memory of nonlocal pseudopotential arrays
  omemproj=real(nprojel,kind=8)*8.d0

  ! Work arrays.
  call memspace_work_arrays_sumrho(lr, mworkrho)
  call memspace_work_arrays_locham(lr, mworkham) !n(m)
  ! pot_ion, rhopot, potxc
  nden=3.d0
  ! In Hamiltonian application: pot + psir + work arrays
  npotham=1.d0+nspinor+real(mworkham * 8 * nspinor, kind=8) / omempot
  ! In sumrho: Rho_p + psir + work arrays
  npotden=1.d0+1.d0+real(mworkrho * 8, kind=8) / omempot
  ! Mixing arrays.
  if (itrpmax /= 1) then
     if (mod(iscf, 10) == 1) narr = 5
     if (mod(iscf, 10) == 2) narr = 2
     if (mod(iscf, 10) == 3) narr = 3
     if (mod(iscf, 10) == 4) narr = 5
     if (mod(iscf, 10) == 5) narr = 10
     if (mod(iscf, 10) == 6) narr = 10
     if (mod(iscf, 10) == 7) narr = 1 + 2 * 7
     nden = nden + narr
  end if

  mem%submat = real(norb,kind=8)**2
  mem%norb = norb
  mem%oneorb = omemwf
  mem%ncomponents = lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
  !takes into account psit
  if(nproc > 1 ) mem%allpsi_mpi=24.d0*real(norbp*nvctrp*nproc,kind=8)
  if(nproc == 1 ) mem%allpsi_mpi=16.d0*real(norbp*nvctrp*nproc,kind=8)
  mem%norbp = norbp
  if(nproc > 1 ) mem%psistorage=8.d0*real(2*idsx+3,kind=8)*real(norbp*nvctrp*nproc,kind=8)
  if(nproc == 1 ) mem%psistorage=8.d0*real(2*idsx+2,kind=8)*real(norbp*nvctrp*nproc,kind=8)
  mem%projarr = omemproj
  mem%grid = omempot
  mem%workarr = real(max(mworkrho,mworkham),kind=8)

  if (nproc > 1) then 
     mem%kernel=19.d0*omemker
     mem%density=mem%psistorage+nden*omemden+npotden*omempot+omemker+omemproj
     mem%psolver=12.d0*omemden+mem%psistorage+omemker+omemproj
     mem%ham=nden*omemden+npotham*omempot+mem%psistorage+omemker+omemproj
  else
     mem%kernel=11.d0*omemker
     mem%density=mem%psistorage+nden*omemden+(npotden-1.d0)*omempot+omemker+omemproj
     mem%psolver=8.d0*omemden+mem%psistorage+omemker+omemproj
     mem%ham=nden*omemden+(npotham-1.d0)*omempot+mem%psistorage+omemker+omemproj
  end if
  !estimation of the memory peak
  mem%peak=max(mem%kernel,mem%density,mem%psolver,mem%ham+mem%submat)

END SUBROUTINE MemoryEstimator

!> old timing routine, should disappear as soon as the f_timing routine is called
subroutine timing(comm,category,action)
  use dictionaries, only: max_field_length,f_err_raise
  use yaml_strings, only: yaml_toa
  use module_types, only: find_category
  use time_profiling
  use wrapper_MPI
  implicit none
  !Variables
  integer, intent(in) :: comm
  character(len=*), intent(in) :: category
  character(len=2), intent(in) :: action  
  !Local variables
  integer :: cat_id,ierr
  character(len=max_field_length) :: cattmp
  external :: gather_timings
  !this is to ensure that timing routines have been properly called
  !call check_initialization()

  !modification of the timing to see if it works
  select case(action)
  case('PR')
     !here iproc is the communicator
     call f_timing_checkpoint(ctr_name=category,mpi_comm=comm,nproc=mpisize(comm),&
          gather_routine=gather_timings)
  case default
     !find category in the old scheme
     call find_category(category,cat_id)

     if (cat_id /= TIMING_UNINITIALIZED) then
        call get_category_name(cat_id,cattmp)
        if (f_err_raise(trim(cattmp)/=trim(category),'Error in category '//&
             trim(yaml_toa(cat_id))//' (name='//trim(category)//' ), found '//&
             trim(cattmp)//' instead',err_name='BIGDFT_RUNTIME_ERROR')) return
        
        call f_timing(cat_id,action)
     end if
  end select

END SUBROUTINE timing

