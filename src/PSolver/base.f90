!> @file
!!    Define types for reals, memory and timing
!! @author 
!!    Copyright (C) 2008-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!>  Modules which contains the low level definitions, as well as some profiling procedures
!!  in order to have a standalone PSolver version.
module module_base 
  !use MPI
  implicit none  
  !buffer to be added at the end of the last dimension of an array to control bounds_check
  integer, parameter :: ndebug=0

  !general precision, density and the wavefunctions types
  integer, parameter :: gp=kind(1.0d0)  !general-type precision
  integer, parameter :: dp=kind(1.0d0)  !density-type precision
  integer, parameter :: wp=kind(1.0d0)  !wavefunction-type precision

  !MPI definitions and datatypes for density and wavefunctions
  include 'mpif.h'
  integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION,mpidtypd=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
  !integer, parameter :: mpidtypw=MPI_REAL,mpidtypd=MPI_REAL !in case of single precision

  !flag for GPU computing, if CUDA libraries are present
  !the parameter flag is removed to allow each processor to modify it following the needs
  !in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  logical :: GPUconv=.false.,GPUblas=.false.

  type :: memstat
     character(len=36) :: routine,array
     integer(kind=8) :: memory,peak
  end type memstat

  ! Save values for memocc.
  type(memstat) :: loc,tot
  integer :: nalloc,ndealloc,iproc

  !interface for the memory allocation control, depends on ndebug
  interface memocc
     module procedure mo_dp1,mo_dp2,mo_dp3,mo_dp4,mo_dp5,mo_dp6,mo_dp7,&
          mo_sp1,mo_sp2,mo_sp3,mo_sp4,mo_sp5,mo_sp6,mo_sp7,&
          mo_i1,mo_i2,mo_i3,mo_i4,mo_i5,mo_i6,mo_i7,&
          mo_l1,mo_l2,mo_l3,mo_l4,mo_l5,mo_l6,mo_l7,&
          mo_c1, mo_cmpdp1,&
          memocc_internal  !central routine to be used for deallocation
  end interface

  contains

    !routine used for deallocations
    subroutine memocc_internal(istat,isize,array,routine)
      implicit none
      character(len=*), intent(in) :: array,routine
      integer, intent(in) :: istat,isize
      call memory_occupation(istat,isize,array,routine) !this routine is in profiling/memory.f90
    END SUBROUTINE memocc_internal

    subroutine dp_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      real(kind=8), dimension(*) :: array
      !local variables
      integer :: i
      real(kind=8), external :: d_nan
      do i=1,npaddim*ndebug
         array(nstart+i)= d_nan() !this function is in profiling/memory.f90
      end do
    END SUBROUTINE dp_padding

    subroutine sp_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      real(kind=4), dimension(*) :: array
      !local variables
      integer :: i
      real(kind=4), external :: r_nan
      do i=1,npaddim*ndebug
         array(nstart+i)= r_nan() !this function is in profiling/memory.f90
      end do
    END SUBROUTINE sp_padding

    subroutine i_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      integer, dimension(*) :: array
      !local variables
      integer :: i
      real(kind=4), external :: r_nan
      do i=1,npaddim*ndebug
         array(nstart+i)= r_nan() !this function is in profiling/timem.f90
      end do
    END SUBROUTINE i_padding

    subroutine l_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      logical, dimension(*) :: array
      !local variables
      integer :: i
      do i=1,npaddim*ndebug
         array(nstart+i)=.false.
      end do
    END SUBROUTINE l_padding

    subroutine c_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      character(len=20), dimension(*) :: array
      !local variables
      integer :: i
      do i=1,npaddim*ndebug
         array(nstart+i)='AAAAAAAAAAAAAAAAAAAA'
      end do
    END SUBROUTINE c_padding

    !beginning of the verbose section
    subroutine mo_dp1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(shape(array))-ndebug
         call dp_padding(1,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_dp1

    subroutine mo_dp2(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:1))
         ndim=product(shape(array))-ndebug*npaddim
         call dp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_dp2

    subroutine mo_dp3(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:2))
         ndim=product(shape(array))-ndebug*npaddim
         call dp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_dp3

    subroutine mo_dp4(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:3))
         ndim=product(shape(array))-ndebug*npaddim
         call dp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_dp4

    subroutine mo_dp5(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:4))
         ndim=product(shape(array))-ndebug*npaddim
         call dp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_dp5

    subroutine mo_dp6(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:5))
         ndim=product(shape(array))-ndebug*npaddim
         call dp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_dp6

    subroutine mo_dp7(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:6))
         ndim=product(shape(array))-ndebug*npaddim
         call dp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_dp7

    subroutine mo_sp1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(shape(array))-ndebug
         call sp_padding(1,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_sp1

    subroutine mo_sp2(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:1))
         ndim=product(shape(array))-ndebug*npaddim
         call sp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_sp2

    subroutine mo_sp3(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:2))
         ndim=product(shape(array))-ndebug*npaddim
         call sp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_sp3

    subroutine mo_sp4(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:3))
         ndim=product(shape(array))-ndebug*npaddim
         call sp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_sp4

    subroutine mo_sp5(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:4))
         ndim=product(shape(array))-ndebug*npaddim
         call sp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_sp5

    subroutine mo_sp6(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:5))
         ndim=product(shape(array))-ndebug*npaddim
         call sp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_sp6

    subroutine mo_sp7(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:6))
         ndim=product(shape(array))-ndebug*npaddim
         call sp_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_sp7

    subroutine mo_i1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(shape(array))-ndebug
         call i_padding(1,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_i1

    subroutine mo_i2(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:1))
         ndim=product(shape(array))-ndebug*npaddim
         call i_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_i2

    subroutine mo_i3(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:2))
         ndim=product(shape(array))-ndebug*npaddim
         call i_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_i3

    subroutine mo_i4(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:3))
         ndim=product(shape(array))-ndebug*npaddim
         call i_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_i4

    subroutine mo_i5(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:4))
         ndim=product(shape(array))-ndebug*npaddim
         call i_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_i5

    subroutine mo_i6(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:5))
         ndim=product(shape(array))-ndebug*npaddim
         call i_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_i6

    subroutine mo_i7(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:6))
         ndim=product(shape(array))-ndebug*npaddim
         call i_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_i7

    subroutine mo_l1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(shape(array))-ndebug
         call l_padding(1,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_l1

    subroutine mo_l2(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:1))
         ndim=product(shape(array))-ndebug*npaddim
         call l_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_l2

    subroutine mo_l3(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:2))
         ndim=product(shape(array))-ndebug*npaddim
         call l_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_l3

    subroutine mo_l4(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:3))
         ndim=product(shape(array))-ndebug*npaddim
         call l_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_l4

    subroutine mo_l5(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:4))
         ndim=product(shape(array))-ndebug*npaddim
         call l_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_l5

    subroutine mo_l6(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:5))
         ndim=product(shape(array))-ndebug*npaddim
         call l_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_l6

    subroutine mo_l7(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=shape(array)
         npaddim=product(iashp(1:6))
         ndim=product(shape(array))-ndebug*npaddim
         call l_padding(npaddim,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_l7

    subroutine mo_c1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      character(len=20), dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(shape(array))-ndebug
         call c_padding(1,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    END SUBROUTINE mo_c1

    subroutine mo_cmpdp1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      complex(kind=8), dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(shape(array))-ndebug
         call cmpdp_padding(1,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_cmpdp1

end module module_base

